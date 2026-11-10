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
from datetime import datetime

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

#Define Caches
POSITIVE_MATCH_CACHE = {}
NEGATIVE_MATCH_CACHE = set()

# Paths (update as needed; can be overridden via CLI)
PDF_DIR = r"D:/data/Max/stress-biomarkers-in-entometabolomics/pdfs"
HMDB_PATH = r"D:/data/Max/stress-biomarkers-in-entometabolomics/hmdb_metabolites/hmdb_metabolites.xml"
LIPDMAPS_PATH = r"D:/data/Max/stress-biomarkers-in-entometabolomics/lipidmaps/lipidmaps.sdf"
TAXON_PATH = r"D:/data/Max/stress-biomarkers-in-entometabolomics/Taxon/Taxon.txt"
LOG_PATH = "skipped_papers.log"
OUTPUT_XLSX = "metabolite_results.xlsx"

# Tunable extraction parameters (can be overridden via CLI)
MAX_NGRAM = 2
TOKEN_MIN_LEN = 4
FUZZY_THRESHOLD = 90
FUZZY_LEN_DIFF_RATIO = 0.2
FUZZY_MIN_LEN = 15
FUZZY_MAX_LEN = 30
FUZZY_K = 2000
FUZZY_M = 30

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

def load_fuzzy_exceptions(path="fuzzy_exceptions.txt"):
    with open(path, encoding="utf-8") as f:
        return {
            normalize(line.strip())
            for line in f
            if line.strip()
        }

EXCLUDED_FUZZY_QUERIES = load_fuzzy_exceptions()

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
            "HMDB/LMP_ID": accession,
            "SMILES": smiles,
            "INCHIKEY": inchikey,
            "CHEMICAL_FORMULA": formula,
            "AVERAGE_MASS": avg_mass,
            "MONO_MASS": mono_mass,
            "HMDB_SUPER_CLASS": hmdb_super_class,
            "HMDB_CLASS": hmdb_class,
            "HMDB_SUB_CLASS": hmdb_sub_class,
            "DATABASE_SOURCE": "HMDB",
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

        # Free memory
        elem.clear()

    return db

def load_lipidmaps(path):
    db = {}

    with open(
        path,
        encoding="utf-8",
        errors="ignore"
    ) as f:
        content = f.read()
    records = content.split("$$$$")
    for record in records:
        record = record.strip()
        if not record:
            continue
        fields = {}
        current = None
        buffer = []
        for line in record.splitlines():
            m = re.match(
                r">\s*<(.+?)>",
                line
            )
            if m:
                if current:
                    fields[current] = (
                        "\n".join(buffer).strip()
                    )
                current = m.group(1)
                buffer = []
            elif current:
                buffer.append(line)
        if current:
            fields[current] = (
                "\n".join(buffer).strip()
            )
        lm_id = fields.get("LM_ID")
        if not lm_id:
            continue
        row = {
            "NAME":
                fields.get("NAME", ""),
            "HMDB/LMP_ID":
                lm_id,
            "SMILES":
                fields.get("SMILES", ""),
            "INCHIKEY":
                fields.get("INCHI_KEY",""),
            "CHEMICAL_FORMULA":
                fields.get("FORMULA",""),
            "AVERAGE_MASS":
                fields.get("EXACT_MASS",""),
            "MONO_MASS":
                fields.get("EXACT_MASS",""),
            "HMDB_SUPER_CLASS":
                "",
            "HMDB_CLASS":
                fields.get("MAIN_CLASS",""),
            "HMDB_SUB_CLASS":
                fields.get("CATEGORY",""),
            "DATABASE_SOURCE":
                "LIPIDMAPS",
            "ABBREVIATION":fields.get("ABBREVIATION","")
        }
        if row["NAME"]:
            db.setdefault(normalize(row["NAME"]),row)
        systematic = fields.get("SYSTEMATIC_NAME", "")
        if systematic:
            db.setdefault(normalize(systematic), row)
        abbr = fields.get("ABBREVIATION","")
        if abbr:
            db.setdefault(normalize(abbr),row)
        synonyms = fields.get("SYNONYMS", "")
        for syn in re.split(r"[;\n|]",synonyms):
            syn = syn.strip()
            if syn:
                db.setdefault(normalize(syn), row)
    return db

def build_combined_db(
    hmdb_db,
    lipid_db):
    combined = hmdb_db.copy()
    added = 0
    for key, value in lipid_db.items():
        if key not in combined:

            combined[key] = value
            added += 1
    return combined

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


COMMON_STRESS_WORDS = set(['mass','produce','fects','tive','induced','cross','facilitated','factors','better','high','higher','strong','weak','under','vation','increased','decreased','tion','its','it','itself','have','has','had','is','are','was','were','be','been','being','do','does','did','will','shall','can','could','may','might','must','should','the', 'a', 'in', 'to', 'of', 'for', 'by', 'with', 'as', 'on', 'at', 'from', 'and', 'or', 'nor', 'but', 'so', 'yet', 'all', 'any', 'during', 'meets', 'multiple', 'post', 'subsequent', 'when', 'adverse', 'affects', 'all', 'although', 'among', 'analysis', 'avoid', 'broadens', 'causes', 'causing', 'confers', 'controlling', 'critical', 'dative', 'different', 'greater', 'hardening', 'ical', 'idative', 'improved', 'induce', 'induces', 'lution', 'mal', 'mediate', 'otic', 'pathways', 'severe', 'showed', 'stronger', 'sublethal', 'term', 'their', 'that', 'upon', 'various', 'lution', 'avoid'])


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


def process_pdf(pdf_path, combined_db, insect_db, genus_to_species_counter, cfg):
    pdf_name = os.path.basename(pdf_path)
    start_time = datetime.now()
    print(f"[{start_time:%d-%m-%Y %H:%M}]" f"Starting PDF: {pdf_name}")
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
        sections.get('discussion', ''),
        sections.get('conclusion', '')
            ]
    relevant_metabolite_text = ' '.join(metabolite_sections)

    tokens = extract_metabolite_tokens(relevant_metabolite_text, token_min_len=cfg['token_min_len'])
    ngrams = generate_ngrams(tokens, max_n=cfg['max_n'])
    norm_ngrams = set(normalize(ng) for ng in ngrams)

    matched = []
    matched_norm_ngrams = set()
    seen_hmdb_ids = set()
    for ng in norm_ngrams:
        if ng in POSITIVE_MATCH_CACHE:
            entry = POSITIVE_MATCH_CACHE[ng]
            hmdb_id = entry["HMDB/LMP_ID"]
            print(f"[POS CACHE HIT] {ng}")

            if hmdb_id not in seen_hmdb_ids:
                seen_hmdb_ids.add(hmdb_id)
                matched.append(entry)
            continue    
        # Already known to never match
        if ng in NEGATIVE_MATCH_CACHE:
            print(f"[NEG CACHE HIT] {ng}")
            continue

        # EXACT MATCH
        if ng in combined_db:

            entry = combined_db[ng]

            name = normalize(entry["NAME"])

            if name in EXCLUDED_METABOLITES:
                continue

            hmdb_id = entry["HMDB/LMP_ID"]

            if hmdb_id in seen_hmdb_ids:
                continue

            seen_hmdb_ids.add(hmdb_id)

            POSITIVE_MATCH_CACHE[ng] = entry
            entry_copy = entry.copy()
            entry_copy["MATCH_TYPE"] = "EXACT"

            matched.append(entry_copy)
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
        candidate_norm = normalize(candidate)
        if candidate_norm in EXCLUDED_FUZZY_QUERIES:
            print(f"[FUZZY EXCLUDED] {candidate}")
            continue
        if fuzzy_match_count >= cfg['fuzzy_m']:
            break  
        if candidate in POSITIVE_MATCH_CACHE:
            entry = POSITIVE_MATCH_CACHE[candidate]

            hmdb_id = entry["HMDB/LMP_ID"]
            if hmdb_id not in seen_hmdb_ids:
                seen_hmdb_ids.add(hmdb_id)
                matched.append(entry)
            print(f"[POS CACHE HIT] {candidate}")

            continue 
          
        candidate_len = len(candidate)
        
        result = process.extractOne(
            candidate,
            METABOLITE_KEYS,
            scorer=fuzz.WRatio,
            score_cutoff=cfg['fuzzy_threshold']
        )

        if not result:
            NEGATIVE_MATCH_CACHE.add(candidate)
            continue 

        best_match, score, _ = result

        hmdb_len = len(best_match)

        len_diff_ratio = (
            abs(candidate_len - hmdb_len)
            / max(candidate_len, hmdb_len)
        )

        if len_diff_ratio > cfg['fuzzy_len_diff_ratio']:
            continue
        if not best_match:
            NEGATIVE_MATCH_CACHE.add(candidate)
            continue

        entry = combined_db[best_match]

        name = normalize(entry["NAME"])

        if name in EXCLUDED_METABOLITES:
            continue
        hmdb_id = entry["HMDB/LMP_ID"]

        if hmdb_id in seen_hmdb_ids:
            continue
        print(
            f"[FUZZY]"
            f"'{candidate}'->"
            f"'{best_match}' "
            f"(score={score})"
        )
        POSITIVE_MATCH_CACHE[candidate] = entry
        
        seen_hmdb_ids.add(hmdb_id)

        print(f"[POS CACHE ADD] {candidate}")
        entry_copy = entry.copy()
        
        entry_copy["MATCH_TYPE"] = "FUZZY"
        entry_copy["FUZZY_QUERY"] = candidate
        entry_copy["FUZZY_MATCH"] = best_match
        
        matched.append(entry_copy)
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
        print("HMDB/LMP ID:", m.get('HMDB/LMP_ID'))
        print("HMDB super_class:", m.get('HMDB_SUPER_CLASS'))
        print("HMDB class:", m.get('HMDB_CLASS'))
        print("HMDB sub_class:", m.get('HMDB_SUB_CLASS'))
        logp = get_logp(m.get('INCHIKEY', ''))
        results.append({
    'Metabolite name': m.get('NAME', ''),
    'InChIKey': m.get('INCHIKEY', ''),
    'HMDB ID/LMP': m.get('HMDB/LMP_ID', ''),
    'CHEMICAL_FORMULA': m.get('CHEMICAL_FORMULA', ''),
    'Average Mass': m.get('AVERAGE_MASS', ''),
    'Mono Mass': m.get('MONO_MASS', ''),
    'LogP': logp,
    'HMDB super_class': m.get('HMDB_SUPER_CLASS', ''),
    'HMDB class': m.get('HMDB_CLASS', ''),
    'HMDB sub_class': m.get('HMDB_SUB_CLASS', ''),
    'Database': m.get('DATABASE_SOURCE',''),
    'Match type': m.get('MATCH_TYPE',''),
    'Fuzzy query': m.get('FUZZY_QUERY',''),
    'Fuzzy match': m.get('FUZZY_MATCH'),
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
    p.add_argument("--lipidmaps-path",default=LIPDMAPS_PATH)

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
    LIPIDMAPS_PATH = args.lipidmaps_path
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
    lipid_db = {}
    if args.lipidmaps_path:
        lipid_db = load_lipidmaps(
            args.lipidmaps_path
        )
    combined_db = build_combined_db(
        hmdb_db,
        lipid_db
    )
    global METABOLITE_KEYS
    METABOLITE_KEYS = list(combined_db.keys())
    
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
            res = process_pdf(pdf, combined_db, insect_db, genus_to_species_counter, cfg)
            if res:
                all_results.extend(res)
        except Exception as e:
            logging.exception(f"Error processing {pdf}: {e}")

    if all_results:
        print("RAW RESULTS", len(all_results))
        merged = {}
        for row in all_results:
            key = (row.get('InChIKey') or row.get('Metabolite name'))
            if key not in merged:
                merged[key] = row.copy()
            else:
                for k, v in row.items():
                    if k in ['stress sources', 'insect species', 'techniques','DOI references']:
                        old = set(merged[key][k].split(', ')) if merged[key][k] else set()
                        new = set(v.split(', ')) if v else set()
                        merged[key][k] = ', '.join(sorted(s.strip() for s in (old | new) if s.strip()))

        df = pd.DataFrame(all_results)
        print("MERGED RESULTS" ,len(df))
        df["n_papers"] = df["DOI references"].apply(
            lambda x: len(
                {
                    doi.strip()
                    for doi in str(x).split(",")
                    if doi.strip()
                }
            )
        )

        per_pdf = {}
        for row in all_results:
            pdf = row.get('DOI references', '').strip()
            if pdf not in per_pdf:
                per_pdf[pdf] = {
                    'metabolites': set(),
                    'techniques': set(),
                    'insect species': set(),
                    'stress sources': set(),
                    'HMDB super_class':set(),
                    'HMDB class':set(),
                    'HMDB sub_class':set(),
                }

            mname = row.get('Metabolite name', '').strip()
            if (
                mname
                and normalize(mname) not in EXCLUDED_METABOLITES
            ):
                per_pdf[pdf]['metabolites'].add(mname)
                super_class = row.get('HMDB super_class','.').strip()
                if super_class:
                    per_pdf[pdf]['HMDB super_class'].add(super_class)
                HMDB_class = row.get('HMDB class','.').strip()
                if HMDB_class:
                    per_pdf[pdf]['HMDB class'].add(HMDB_class)
                sub_class = row.get('HMDB sub_class','.').strip()
                if sub_class:
                    per_pdf[pdf]['HMDB sub_class'].add(sub_class)
                
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

        per_pdf_rows = []
        for pdf, vals in per_pdf.items():
            per_pdf_rows.append({
                'DOI references':pdf,
                'n_metabolites': len(vals['metabolites']),
                'metabolites': '; '.join(sorted(vals['metabolites'])),
                'HMDB super_class': '; '.join(sorted(vals['HMDB super_class'])),
                'HMDB class': '; '.join(sorted(vals['HMDB class'])),
                'HMDB sub_class': '; '.join(sorted(vals['HMDB sub_class'])),
                'techniques': ', '.join(sorted(vals['techniques'])),
                'insect species': ', '.join(sorted(vals['insect species'])),
                'stress sources': ', '.join(sorted(vals['stress sources'])),
                })

        df_pdf_summary = pd.DataFrame(per_pdf_rows).sort_values('DOI references')

        with pd.ExcelWriter(OUTPUT_XLSX, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='metabolite_hits', index=False)
            df_pdf_summary.to_excel(writer, sheet_name='per_pdf_summary', index=False)
        print(f"Results saved to {OUTPUT_XLSX}")
    else:
        print("No metabolites extracted from any PDFs.")
   
if __name__ == '__main__':
    main()
