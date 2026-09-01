import os
import re
import csv
import sys
import logging
import unicodedata
import urllib.parse
import argparse
from collections import defaultdict, Counter
from functools import lru_cache

import pandas as pd
from pdfminer.high_level import extract_text
import requests
from rapidfuzz import fuzz, process
import xml.etree.ElementTree as ET
import pickle

csv.field_size_limit(sys.maxsize)

# Constants
HMDB_DELIMITER_SAMPLE_SIZE = 8192  # Bytes to read for delimiter detection
HMDB_MIN_SAMPLE_SIZE = 100  # Minimum bytes required for reliable delimiter detection
HMDB_DELIMITER = ','  # Force comma delimiter for HMDB metabolite files

# Paths (update as needed; can be overridden via CLI)
PDF_DIR = r"D:/data/Max/stress-biomarkers-in-entometabolomics/pdfs_crossval"
HMDB_PATH = r"D:/data/Max/stress-biomarkers-in-entometabolomics/hmdb_metabolites/hmdb_metabolites.xml"
TAXON_PATH = r"D:/data/Max/stress-biomarkers-in-entometabolomics/Taxon/Taxon.txt"
LOG_PATH = "skipped_papers.log"
OUTPUT_XLSX = "metabolite_results_fuzzy.xlsx"

# Tunable extraction parameters (can be overridden via CLI)
MAX_NGRAM = 2
TOKEN_MIN_LEN = 3
FUZZY_THRESHOLD = 70
FUZZY_LEN_DIFF_RATIO = 0.5
FUZZY_MIN_LEN = 2
FUZZY_MAX_LEN = 80
FUZZY_K = 12000
FUZZY_M = 400

# List or set of PDF filenames for which to enable debug output
# Add PDF filenames (e.g. {'file1.pdf', 'file2.pdf'}) to enable debug for specific files
# Leave as empty set for no debug output

debug_pdfs = set()

# Setup logging
logging.basicConfig(filename=LOG_PATH, level=logging.INFO, format='%(asctime)s %(message)s')

# Pre-compiled regex patterns for better performance
# Updated to preserve lipid notation characters: keep : / ( ) . , and hyphens inside tokens
# Only remove em-dashes, en-dashes, semicolons, and collapse whitespace
_normalize_re = re.compile(r'[–—;]|\s+')
_word_re = re.compile(r'\b\w+\b')
_hyphen_break_re = re.compile(r'-\s*\n\s*')
_doi_re = re.compile(r"10\.\d{4,9}/[-._;()/:A-Z0-9]+", re.IGNORECASE)
_element_re = re.compile(r'[A-Z][a-z]?')
_stress_re = re.compile(r'\b(?:\w+\s)?(?:stress|tolerance)\b', re.IGNORECASE)
_insect_species_re = re.compile(r'\b[A-Z][a-z]{2,}\s*[ \n\r\t]+\s*[a-z]{3,}\b')
_abbrev_species_re = re.compile(r'\b([A-Z])\.\s*([a-z]{3,})\b')
_whitespace_re = re.compile(r'\s+')
_non_alpha_space_re = re.compile(r'[^a-z ]')
_parentheses_re = re.compile(r'\s*\(.*?\)\s*')
_synonym_separator_re = re.compile(r'[|;,\t]')
# Metabolite token regex: captures tokens with digits, hyphens, primes, and chemical punctuation (/:(),.)
_metabolite_token_re = re.compile(r'\b[\w][\w\-′″‴⁗\'\"/:\(\),\.]*[\w]\b|\b\w\b')

# Translation table for Greek letters and primes
_greek_translation = str.maketrans({
    'α': 'alpha', 'β': 'beta', 'γ': 'gamma', 'δ': 'delta', 'ε': 'epsilon',
    'ζ': 'zeta', 'η': 'eta', 'θ': 'theta', 'ι': 'iota', 'κ': 'kappa',
    'λ': 'lambda', 'μ': 'mu', 'ν': 'nu', 'ξ': 'xi', 'ο': 'omicron',
    'π': 'pi', 'ρ': 'rho', 'σ': 'sigma', 'ς': 'sigma', 'τ': 'tau',
    'υ': 'upsilon', 'φ': 'phi', 'χ': 'chi', 'ψ': 'psi', 'ω': 'omega',
    'Α': 'alpha', 'Β': 'beta', 'Γ': 'gamma', 'Δ': 'delta', 'Ε': 'epsilon',
    'Ζ': 'zeta', 'Η': 'eta', 'Θ': 'theta', 'Ι': 'iota', 'Κ': 'kappa',
    'Λ': 'lambda', 'Μ': 'mu', 'Ν': 'nu', 'Ξ': 'xi', 'Ο': 'omicron',
    'Π': 'pi', 'Ρ': 'rho', 'Σ': 'sigma', 'Τ': 'tau', 'Υ': 'upsilon',
    'Φ': 'phi', 'Χ': 'chi', 'Ψ': 'psi', 'Ω': 'omega',
    '′': '', '″': '', '‴': '', '⁗': '', "'": '', '"': ''
})

# Technique/platform patterns
TECHNIQUE_PATTERNS = {
    "untargeted metabolomics": re.compile(r"\b(?:un[\-\s]?targeted|non[\-\s]?targeted)\s+metabolomics\b", re.IGNORECASE),
    "targeted metabolomics": re.compile(r"\btargeted\s+metabolomics\b", re.IGNORECASE),
    "LC-MS": re.compile(r"\b(?:U?H?PLC|LC)[\-\–\—\s]?MS(?:/MS)?\b|\b(?:LC|U?H?PLC)[\-\–\—\s]?(?:HRMS|high[\-\–\—\s]?resolution mass spectrometr(?:y|ic))\b|\bliquid chromatography[\-\–\—\s]*(?:mass spectrometry|mass spectrometr(?:y|ic))\b", re.IGNORECASE),
    "GC-MS": re.compile(r"\bGC[\-\–\—\s]?MS(?:/MS)?\b|\bGC[\-\–\—\s]?(?:HRMS|high[\-\–\—\s]?resolution mass spectrometr(?:y|ic))\b|\bgas chromatography[\-\–\—\s]*(?:mass spectrometry|mass spectrometr(?:y|ic))\b", re.IGNORECASE),
    "NMR": re.compile(r"\bNMR\b|\bnuclear magnetic resonance\b", re.IGNORECASE),
}


@lru_cache(maxsize=4096)
def normalize(text):
    text = text.lower().translate(_greek_translation)
    text = _normalize_re.sub(' ', text)
    return ' '.join(text.split())

def load_exceptions(path="exceptions.txt"):
    with open(path, encoding="utf-8") as f:
        return {
            normalize(line.strip())
            for line in f
            if line.strip()
        }

EXCLUDED_METABOLITES = load_exceptions()
def extract_words(text):
    text = _hyphen_break_re.sub('', text)
    words = _word_re.findall(text)
    return [w for w in words if len(w) >= 3]
    
def _safe_text(elem, path, default=''):
    node = elem.find(path)
    if node is None or node.text is None:
        return default
    return node.text.strip()

def extract_metabolite_tokens(text, token_min_len=TOKEN_MIN_LEN):
    text = _hyphen_break_re.sub('', text)
    tokens = _metabolite_token_re.findall(text)
    return [t for t in tokens if len(t) >= token_min_len]


def generate_ngrams(tokens, max_n=MAX_NGRAM):
    leading_stop_tokens = {"the", "and", "of", "in", "for", "with", "by", "on"}

    ngrams = []
    for i in range(len(tokens)):
        first_token_lower = tokens[i].lower()
        if first_token_lower in leading_stop_tokens:
            continue

        for n in range(1, min(max_n + 1, len(tokens) - i + 1)):
            ngram = ' '.join(tokens[i:i+n])
            ngrams.append(ngram)
    return ngrams


def extract_doi(text):
    return _doi_re.findall(text)


def is_single_element(formula):
    return len(_element_re.findall(formula)) == 1


def load_hmdb(path):
    """
    Load HMDB metabolites from XML.
    Returns:
        dict: normalized name/synonym -> metabolite row dict
    """

    db = {}

    HMDB_NS = "{http://www.hmdb.ca}"

    count = 0

    for event, elem in ET.iterparse(path, events=("end",)):

        # HMDB uses namespaces
        if elem.tag != f"{HMDB_NS}metabolite":
            continue

        count += 1

        def get_text(parent, tag):
            node = parent.find(f"{HMDB_NS}{tag}")
            return node.text.strip() if node is not None and node.text else ""

        name = get_text(elem, "name")
        accession = get_text(elem, "accession")
        smiles = get_text(elem, "smiles")
        inchikey = get_text(elem, "inchikey")
        formula = get_text(elem, "chemical_formula")
        avg_mass = get_text(elem, "average_molecular_weight")
        mono_mass = get_text(elem, "monisotopic_molecular_weight")

        taxonomy = elem.find(f"{HMDB_NS}taxonomy")

        hmdb_super_class = ""
        hmdb_class = ""
        hmdb_sub_class = ""

        if taxonomy is not None:
            hmdb_super_class = get_text(taxonomy, "super_class")
            hmdb_class = get_text(taxonomy, "class")
            hmdb_sub_class = get_text(taxonomy, "sub_class")

        row = {
            "NAME": name,
            "HMDB_ID": accession,
            "SMILES": smiles,
            "INCHIKEY": inchikey,
            "CHEMICAL_FORMULA": formula,
            "AVERAGE_MASS": avg_mass,
            "MONO_MASS": mono_mass,
            "HMDB_SUPER_CLASS": hmdb_super_class,
            "HMDB_CLASS": hmdb_class,
            "HMDB_SUB_CLASS": hmdb_sub_class,
        }

        # Primary name
        if name:
            db.setdefault(normalize(name), row)

        # Synonyms
        synonyms = elem.find(f"{HMDB_NS}synonyms")
        if synonyms is not None:

            for syn in synonyms.findall(f"{HMDB_NS}synonym"):

                if syn.text:
                    syn_key = normalize(syn.text.strip())

                    if syn_key:
                        db.setdefault(syn_key, row)

        # Debug first few metabolites only
        if count <= 5:
            print(
                f"[HMDB LOAD] "
                f"{accession} | "
                f"{name} | "
                f"subclass={hmdb_sub_class}"
            )

        # Free memory
        elem.clear()

    return db

def normalize_text(text):
    return unicodedata.normalize("NFKD", text)


@lru_cache(maxsize=4096)
def normalize_insect_name(name):
    name = normalize_text(name)
    name = _whitespace_re.sub(' ', name)
    name = name.strip().lower()
    name = _parentheses_re.sub('', name)
    name = _non_alpha_space_re.sub('', name)
    name = ' '.join(name.split()[:2])
    return name


def load_taxon(path):
    insects = set()
    with open(path, encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row.get('class', '').strip().lower() == 'insecta':
                name = row.get('scientificName', '')
                norm_name = normalize_insect_name(name)
                if norm_name:
                    insects.add(norm_name)
    return insects


_requests_cache = {}


@lru_cache(maxsize=2048)
def get_logp(inchikey):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/XLogP/JSON'
    if url in _requests_cache:
        return _requests_cache[url]
    try:
        r = requests.get(url, timeout=10)
        if r.ok:
            props = r.json()['PropertyTable']['Properties'][0]
            value = props.get('XLogP', None)
            _requests_cache[url] = value
            return value
    except Exception:
        pass
    return None


@lru_cache(maxsize=2048)
def get_npclassifier_superclass(smiles, inchikey):
    result = None
    if smiles:
        url = f"https://npclassifier.ucsd.edu/classify?smiles={urllib.parse.quote(smiles)}"
        if url in _requests_cache:
            return _requests_cache[url]
        try:
            headers = {
                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/125.0.0.0 Safari/537.36"
            }
            r = requests.get(url, headers=headers, timeout=15)
            r.raise_for_status()
            data = r.json()
            if "superclass" in data:
                result = data["superclass"][0] if data["superclass"] else ""
            elif "superclass_results" in data:
                result = data["superclass_results"][0] if data["superclass_results"] else ""
            _requests_cache[url] = result
        except Exception:
            pass
    return result


def extract_sections(text):
    known_headers = {
        'abstract': 'abstract',
        'introduction': 'introduction',
        'materials and methods': 'materialsandmethods',
        'material and methods': 'materialsandmethods',
        'methods': 'materialsandmethods',
        'results': 'results',
        'keywords': 'keywords',
        'discussion': 'discussion',
        'conclusion': 'conclusion',
        'references': 'references',
    }

    merged_variants = ['results and discussion', 'results & discussion', 'results discussion']

    sections = {}
    current_header = 'preamble'
    buffer_lines = []

    def flush():
        if buffer_lines:
            sections[current_header] = '\n'.join(buffer_lines).strip()
            buffer_lines.clear()

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            buffer_lines.append(raw_line)
            continue
        norm = re.sub(r'[^a-z& ]', '', line.lower())
        norm = re.sub(r'\s+', ' ', norm).strip()
        if norm in merged_variants:
            flush()
            current_header = 'results'
            continue
        if norm in known_headers:
            flush()
            current_header = known_headers[norm]
            continue
        if line.isupper() and len(line) >= 3 and len(line.split()) <= 6:
            flush()
            current_header = re.sub(r'\s+', '', line.lower())
            continue
        buffer_lines.append(raw_line)
    flush()

    if 'results' not in sections and 'discussion' in sections:
        sections['results'] = sections['discussion']

    return sections


def extract_insect_species(section_text, insect_db, debug_pdf=None):
    genus_species_in_db = set()
    genus_species_matches = set()

    for match in re.findall(r'\b[A-Z][a-z]{2,}\s*[ \n\r\t]+\s*[a-z]{3,}\b', section_text):
        norm_match = normalize_insect_name(normalize_text(match))
        if norm_match in insect_db:
            genus_species_matches.add(match)
            genus_species_in_db.add(norm_match.split()[0])

    abbrev_matches = re.findall(r'\b([A-Z])\.\s*([a-z]{3,})\b', section_text)
    genus_map = {}
    for name in insect_db:
        genus = name.split()[0] if ' ' in name else name
        if genus:
            genus_map.setdefault(genus[0], set()).add(genus)

    abbrev_matches_found = set()
    for initial, species in abbrev_matches:
        initial = initial.lower()
        for genus in genus_map.get(initial, []):
            norm_abbrev = f"{genus} {species}"
            if norm_abbrev in insect_db:
                abbrev_matches_found.add(f"{genus} {species}")
                genus_species_in_db.add(genus)

    if genus_species_matches:
        return list(genus_species_matches)
    if abbrev_matches_found:
        return list(abbrev_matches_found)
    return []


COMMON_STRESS_WORDS = set(['eukaryotic','certain','against','different','time','this','mass','produce','induced','cross'])


def extract_stress_sources(text):
    matches = _stress_re.findall(text)
    sources = set()
    for m in matches:
        cleaned = m.lower().replace('tolerance', '').replace('stress', '')
        cleaned = cleaned.strip().replace('-', ' ').replace('_', ' ')
        cleaned = ' '.join(cleaned.split())
        if cleaned:
            sources.add(cleaned)
    filtered_sources = [s for s in sources if s not in COMMON_STRESS_WORDS]
    return sorted(filtered_sources)


def extract_techniques(text):
    found = []
    for label, pattern in TECHNIQUE_PATTERNS.items():
        if pattern.search(text):
            found.append(label)
    return sorted(set(found))


@lru_cache(maxsize=128)
def _cached_extract_text(pdf_path):
    return extract_text(pdf_path)


def process_pdf(pdf_path, hmdb_db, insect_db, genus_to_species_counter, cfg):
    pdf_name = os.path.basename(pdf_path)
    text = _cached_extract_text(pdf_path) if pdf_name in debug_pdfs else extract_text(pdf_path)
    sections = extract_sections(text)

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    title_block = ' '.join(lines[:10])

    relevant_insect_text = ' '.join([
        title_block,
        sections.get('abstract', ''),
        sections.get('keywords', ''),
        sections.get('materialsandmethods', '')
    ])

    insect_species = extract_insect_species(relevant_insect_text, insect_db)
    normalized_insect_species = []
    norm_species = [normalize_insect_name(s) for s in insect_species]

    if norm_species and all(' ' not in s for s in norm_species):
        for genus in norm_species:
            if genus in genus_to_species_counter and genus_to_species_counter[genus]:
                most_common = genus_to_species_counter[genus].most_common(1)
                normalized_insect_species.extend([f"{genus} {species}" for species, _ in most_common])
    else:
        normalized_insect_species = sorted(set(norm_species))

    if not normalized_insect_species:
        return None

    metabolite_sections = [
        title_block,
        sections.get('abstract', ''),
        sections.get('results', ''),
        sections.get('keywords', ''),
        sections.get('introduction', ''),
        sections.get('discussion', ''),
        sections.get('conclusion', ''),
        sections.get('references', '')
    ]
    relevant_metabolite_text = ' '.join(metabolite_sections)

    tokens = extract_metabolite_tokens(relevant_metabolite_text, token_min_len=cfg['token_min_len'])
    ngrams = generate_ngrams(tokens, max_n=cfg['max_n'])
    norm_ngrams = set(normalize(ng) for ng in ngrams)

    matched = []
    matched_norm_ngrams = set()
    seen_hmdb_ids = set()
    for ng in norm_ngrams:
        # EXACT MATCH
        if ng in hmdb_db:

            entry = hmdb_db[ng]

            name = normalize(entry["NAME"])

            if (
                name in EXCLUDED_METABOLITES
                ):
                continue

            hmdb_id = entry["HMDB_ID"]

            if hmdb_id in seen_hmdb_ids:
                continue

            seen_hmdb_ids.add(hmdb_id)

            matched.append(entry)
            matched_norm_ngrams.add(ng)

            continue
        
        # FUZZY FALLBACK
    non_matched_ngrams = [
        ng
        for ng in norm_ngrams
        if ng not in matched_norm_ngrams
    ]
    fuzzy_candidates = [
        ng
        for ng in non_matched_ngrams
        if(
            cfg['fuzzy_min_len']
            <= len(ng)
            <=cfg['fuzzy_max_len']
        )
    ]
    fuzzy_candidates = sorted(
        set(fuzzy_candidates),
        key=len
    )[:cfg['fuzzy_k']]

    fuzzy_match_count = 0
    
    for candidate in fuzzy_candidates:
        if fuzzy_match_count >= cfg['fuzzy_m']:
            break
        
        candidate_len = len(candidate)
        
        result = process.extractOne(
            candidate,
            HMDB_KEYS,
            scorer=fuzz.token_ratio,
            score_cutoff=cfg['fuzzy_threshold']
        )

        if not result:
            continue

        best_match, score, _ = result

        hmdb_len = len(best_match)

        len_diff_ratio = (
            abs(candidate_len - hmdb_len)
            / max(candidate_len, hmdb_len)
        )

        if len_diff_ratio > cfg['fuzzy_len_diff_ratio']:
            continue
        entry = hmdb_db[best_match]

        name = normalize(entry["NAME"])

        if name in EXCLUDED_METABOLITES:
            continue
        hmdb_id = entry["HMDB_ID"]

        if hmdb_id in seen_hmdb_ids:
            continue
        print(
            f"[FUZZY]"
            f"'{candidate}'->"
            f"'{best_match}' "
            f"(score={score})"
        )

        seen_hmdb_ids.add(hmdb_id)

        matched.append(entry)
        matched_norm_ngrams.add(best_match)

        fuzzy_match_count += 1
    
    matched = [
        m 
        for m in matched 
        if not is_single_element(
            m.get('CHEMICAL_FORMULA', '')
        )
    ]

    if not matched:
        return None

    relevant_stress_text = ' '.join([
        title_block,
        sections.get('abstract', ''),
        sections.get('results', ''),
        sections.get('keywords', ''),
        sections.get('introduction', ''),
        sections.get('discussion', ''),
        sections.get('conclusion', ''),
        sections.get('references', '')
    ])
    stress_sources = extract_stress_sources(relevant_stress_text)
    techniques = extract_techniques(sections.get('materialsandmethods', ''))

    dois = list(dict.fromkeys(extract_doi(text)))
    doi_reference = dois[0] if dois else ''

    results = []
    for m in matched:
        print("\nMATCH FOUND")
        print("Name:", m.get('NAME'))
        print("HMDB ID:", m.get('HMDB_ID'))
        print("HMDB super_class:", m.get('HMDB_SUPER_CLASS'))
        print("HMDB class:", m.get('HMDB_CLASS'))
        print("HMDB sub_class:", m.get('HMDB_SUB_CLASS'))
        logp = get_logp(m.get('INCHIKEY', ''))
        superclass = get_npclassifier_superclass(m.get('SMILES', ''), m.get('INCHIKEY', ''))
        results.append({
    'source_pdf': pdf_name,
    'Metabolite name': m.get('NAME', ''),
    'InChIKey': m.get('INCHIKEY', ''),
    'HMDB ID': m.get('HMDB_ID', ''),
    'CHEMICAL_FORMULA': m.get('CHEMICAL_FORMULA', ''),
    'Average Mass': m.get('AVERAGE_MASS', ''),
    'Mono Mass': m.get('MONO_MASS', ''),
    'LogP': logp,
    'superclass': superclass,  # from NPClassifier API (keep if useful)
    'HMDB super_class': m.get('HMDB_SUPER_CLASS', ''),
    'HMDB class': m.get('HMDB_CLASS', ''),
    'HMDB sub_class': m.get('HMDB_SUB_CLASS', ''),  # <-- requested
    'stress sources': ', '.join(sorted(set(s for s in stress_sources if s))),
    'insect species': ', '.join(normalized_insect_species),
    'DOI references': doi_reference,
    'techniques': ', '.join(techniques),
})
    return results


def parse_args():
    p = argparse.ArgumentParser(description='Extract stress biomarkers from insect papers')
    p.add_argument('--pdf-dir', default=PDF_DIR)
    p.add_argument('--hmdb-path', default=HMDB_PATH)
    p.add_argument('--taxon-path', default=TAXON_PATH)
    p.add_argument('--log-path', default=LOG_PATH)
    p.add_argument('--output-xlsx', default=OUTPUT_XLSX)

    p.add_argument('--max-n', type=int, default=MAX_NGRAM)
    p.add_argument('--token-min-len', type=int, default=TOKEN_MIN_LEN)
    p.add_argument('--fuzzy-threshold', type=int, default=FUZZY_THRESHOLD)
    p.add_argument('--fuzzy-len-diff-ratio', type=float, default=FUZZY_LEN_DIFF_RATIO)
    p.add_argument('--fuzzy-min-len', type=int, default=FUZZY_MIN_LEN)
    p.add_argument('--fuzzy-max-len', type=int, default=FUZZY_MAX_LEN)
    p.add_argument('--fuzzy-k', type=int, default=FUZZY_K)
    p.add_argument('--fuzzy-m', type=int, default=FUZZY_M)

    return p.parse_args()


def main():
    args = parse_args()

    global PDF_DIR, HMDB_PATH, TAXON_PATH, LOG_PATH, OUTPUT_XLSX
    PDF_DIR = args.pdf_dir
    HMDB_PATH = args.hmdb_path
    TAXON_PATH = args.taxon_path
    LOG_PATH = args.log_path
    OUTPUT_XLSX = args.output_xlsx

    cfg = {
        'max_n': args.max_n,
        'token_min_len': args.token_min_len,
        
        'fuzzy_threshold': args.fuzzy_threshold,
        'fuzzy_len_diff_ratio': args.fuzzy_len_diff_ratio,
        'fuzzy_min_len': args.fuzzy_min_len,
        'fuzzy_max_len': args.fuzzy_max_len,
        
        'fuzzy_k': args.fuzzy_k,
        'fuzzy_m': args.fuzzy_m,
    }

    hmdb_db = load_hmdb(HMDB_PATH)
    global HMDB_KEYS
    HMDB_KEYS = list(hmdb_db.keys())

    print(f"HMDB lookup keys: {len(HMDB_KEYS):,}")
    print("HMDB entries:", len(hmdb_db))
    print(hmdb_db.get(normalize("1-MHis")))
    insect_db = load_taxon(TAXON_PATH)
    pdf_files = [os.path.join(PDF_DIR, f) for f in os.listdir(PDF_DIR) if f.lower().endswith('.pdf')]

    genus_to_species_counter = defaultdict(Counter)
    for pdf in pdf_files:
        text = extract_text(pdf)
        sections = extract_sections(text)
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        title_block = ' '.join(lines[:10])
        relevant_insect_text = ' '.join([
            title_block,
            sections.get('abstract', ''),
            sections.get('keywords', ''),
            sections.get('materialsandmethods', '')
        ])
        for match in re.findall(r'\b[A-Z][a-z]{2,}\s*[ \n\r\t]+\s*[a-z]{3,}\b', relevant_insect_text):
            norm_match = normalize_insect_name(normalize_text(match))
            if norm_match in insect_db and ' ' in norm_match:
                genus, species = norm_match.split()[:2]
                genus_to_species_counter[genus][species] += 1

    all_results = []
    for pdf in pdf_files:
        try:
            res = process_pdf(pdf, hmdb_db, insect_db, genus_to_species_counter, cfg)
            if res:
                all_results.extend(res)
        except Exception as e:
            logging.exception(f"Error processing {pdf}: {e}")

    if all_results:
        merged = {}
        for row in all_results:
            key = (row.get('source_pdf', ''), row.get('InChIKey') or row.get('Metabolite name'))
            if key not in merged:
                merged[key] = row.copy()
            else:
                for k, v in row.items():
                    if k in ['stress sources', 'insect species', 'DOI references', 'techniques']:
                        old = set(merged[key][k].split(', ')) if merged[key][k] else set()
                        new = set(v.split(', ')) if v else set()
                        merged[key][k] = ', '.join(sorted(s.strip() for s in (old | new) if s.strip()))

        df = pd.DataFrame(merged.values())

        # Summary per PDF: metabolites + techniques
        per_pdf = {}
        for row in merged.values():
            pdf = row.get('source_pdf', '')
            if pdf not in per_pdf:
                per_pdf[pdf] = {
                    'source_pdf': pdf,
                    'metabolites': set(),
                    'techniques': set(),
                    'insect species': set(),
                    'stress sources': set(),
                    'DOI references': set(),
                }

            mname = row.get('Metabolite name', '').strip()
            if mname:
                per_pdf[pdf]['metabolites'].add(mname)

            for t in (row.get('techniques', '') or '').split(','):
                t = t.strip()
                if t:
                    per_pdf[pdf]['techniques'].add(t)

            for s in (row.get('insect species', '') or '').split(','):
                s = s.strip()
                if s:
                    per_pdf[pdf]['insect species'].add(s)

            for s in (row.get('stress sources', '') or '').split(','):
                s = s.strip()
                if s:
                    per_pdf[pdf]['stress sources'].add(s)

            for d in (row.get('DOI references', '') or '').split(','):
                d = d.strip()
                if d:
                    per_pdf[pdf]['DOI references'].add(d)

        per_pdf_rows = []
        for pdf, vals in per_pdf.items():
            per_pdf_rows.append({
                'source_pdf': pdf,
                'n_metabolites': len(vals['metabolites']),
                'metabolites': '; '.join(sorted(vals['metabolites'])),
                'techniques': ', '.join(sorted(vals['techniques'])),
                'insect species': ', '.join(sorted(vals['insect species'])),
                'stress sources': ', '.join(sorted(vals['stress sources'])),
                'DOI references': ', '.join(sorted(vals['DOI references'])),
            })

        df_pdf_summary = pd.DataFrame(per_pdf_rows).sort_values('source_pdf')

        with pd.ExcelWriter(OUTPUT_XLSX, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='metabolite_hits', index=False)
            df_pdf_summary.to_excel(writer, sheet_name='per_pdf_summary', index=False)
        print(f"Results saved to {OUTPUT_XLSX}")
    else:
        print("No metabolites extracted from any PDFs.")


if __name__ == '__main__':
    main()
