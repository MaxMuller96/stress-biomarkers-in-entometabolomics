import os
import re
import csv
import sys
import logging
import unicodedata
import urllib.parse
from collections import defaultdict, Counter
from functools import lru_cache

import pandas as pd
from pdfminer.high_level import extract_text
import requests
from rapidfuzz import fuzz

csv.field_size_limit(sys.maxsize)

# Constants
HMDB_DELIMITER_SAMPLE_SIZE = 8192  # Bytes to read for delimiter detection
HMDB_MIN_SAMPLE_SIZE = 100  # Minimum bytes required for reliable delimiter detection
HMDB_DELIMITER = ','  # Force comma delimiter for HMDB metabolite files

# Paths (update as needed)
PDF_DIR = r".../pdfs"
HMDB_PATH = r".../hmdb_metabolites/metabolites-2025-06-20.txt"
TAXON_PATH = r".../Taxon/Taxon.txt"
LOG_PATH = "skipped_papers.log"
OUTPUT_XLSX = "metabolite_results.xlsx"

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
# to preserve complex metabolite names like "N-(2-hydroxyethyl)-palmitamide" or "3,4-dihydroxy-L-phenylalanine"
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

@lru_cache(maxsize=4096)
def normalize(text):
    # Normalize Greek letters and primes to ASCII equivalents using translation table
    text = text.lower().translate(_greek_translation)
    # Replace em-dashes, en-dashes, semicolons with space, then collapse whitespace
    # This preserves lipid notation characters: : / ( ) . , and hyphens within tokens
    text = _normalize_re.sub(' ', text)
    # Remove leading/trailing whitespace and collapse internal whitespace
    return ' '.join(text.split())

def extract_words(text):
    # Remove hyphens at line breaks, normalize, and extract words
    text = _hyphen_break_re.sub('', text)
    words = _word_re.findall(text)
    return [w for w in words if len(w) >= 3]

def extract_metabolite_tokens(text):
    """Extract tokens for metabolite matching, including digits, hyphens, and primes."""
    # Remove hyphens at line breaks
    text = _hyphen_break_re.sub('', text)
    # Use pre-compiled pattern for token extraction
    tokens = _metabolite_token_re.findall(text)
    return [t for t in tokens if len(t) >= 2]

def generate_ngrams(tokens, max_n=8):
    """Generate n-grams from tokens (1-gram to max_n-gram).
    
    Filters applied for optimization:
    - Skip n-grams whose normalized string length exceeds 64 characters
    - Skip n-grams that start with leading stop tokens
    """
    # Define leading stop tokens to skip
    leading_stop_tokens = {"the", "and", "of", "in", "for", "with", "by", "on"}
    
    ngrams = []
    for i in range(len(tokens)):
        for n in range(1, min(max_n + 1, len(tokens) - i + 1)):
            ngram = ' '.join(tokens[i:i+n])
            
            # Skip n-grams that start with leading stop tokens (case-insensitive)
            first_token = tokens[i].lower()
            if first_token in leading_stop_tokens:
                continue
            
            # Skip n-grams whose normalized length exceeds 64 characters
            normalized_ngram = normalize(ngram)
            if len(normalized_ngram) > 64:
                continue
            
            ngrams.append(ngram)
    return ngrams

def extract_doi(text):
    """Extracts all DOI strings from the given text."""
    return _doi_re.findall(text)

def is_single_element(formula):
    # Returns True if the formula contains only one type of element (e.g., 'H2', 'O2')
    return len(_element_re.findall(formula)) == 1


def load_hmdb(path):
    db = {}
    synonym_fields = ['SYNONYMS', 'SYNONYM', 'TRADITIONAL_NAME', 'IUPAC_NAME', 'COMMON_NAME']
    row_count = 0
    
    # Force comma delimiter per provided header
    delimiter = HMDB_DELIMITER
    logging.info(f"HMDB: Using forced delimiter: {repr(delimiter)}")
    
    with open(path, encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            row_count += 1
            # Add primary NAME
            name_norm = normalize(row['NAME'])
            db[name_norm] = row
            
            # Add synonyms and alternative names with multiple separators: |;,\t
            for field in synonym_fields:
                if field in row and row[field]:
                    # Split by multiple delimiters using _synonym_separator_re pattern: [|;,\t]
                    synonyms = _synonym_separator_re.split(row[field])
                    for syn in synonyms:
                        syn = syn.strip()
                        if syn:
                            syn_norm = normalize(syn)
                            if syn_norm and syn_norm not in db:
                                db[syn_norm] = row
    
    unique_name_count = len(db)
    logging.info(f"HMDB: Loaded {row_count} rows, {unique_name_count} unique normalized names")
    return db

def normalize_text(text):
    return unicodedata.normalize("NFKD", text)

@lru_cache(maxsize=4096)
def normalize_insect_name(name):
    name = normalize_text(name)
    name = _whitespace_re.sub(' ', name)  # Collapse all whitespace to a single space
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

from functools import lru_cache

import requests
from functools import lru_cache

# Simple in-memory cache for requests responses
_requests_cache = {}

@lru_cache(maxsize=2048)
def get_logp(inchikey):
    # PubChem REST API
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

import urllib.parse

@lru_cache(maxsize=2048)
def get_npclassifier_superclass(smiles, inchikey):
    # NPClassifier API (GET request with fallback)
    result = None
    if smiles:
        url = f"https://npclassifier.ucsd.edu/classify?smiles={urllib.parse.quote(smiles)}"
        if url in _requests_cache:
            return _requests_cache[url]
        try:
            headers = {
                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/125.0.0.0 Safari/537.36"
            }
            logging.debug(f"NPClassifier GET request: {url}")
            r = requests.get(url, headers=headers, timeout=15)
            r.raise_for_status()
            data = r.json()
            logging.debug(f"NPClassifier GET response: {data}")
            # Handle NPClassifier response format
            if "superclass" in data:
                result = data["superclass"][0] if data["superclass"] else ""
            elif "superclass_results" in data:
                result = data["superclass_results"][0] if data["superclass_results"] else ""
            _requests_cache[url] = result
        except Exception as e:
            logging.exception(f"NPClassifier error: {str(e)}")
    return result

def spaced_re(s):
    # Generate a regex to match section headers with arbitrary whitespace and case-insensitivity
    return ''.join(f'{c}[\s]*' for c in s)

def extract_sections(text):
    """Split the full article text into coarse sections.

    Detection strategy:
    1. Iterate line-by-line; if a line (after stripping) matches one of the known
       headers (case-insensitive, ignoring spaces and punctuation) we start a new
       section.
    2. Lines that are in ALL CAPS and longer than 3 chars are also considered
       potential generic headers (used by some journals).
    3. Merged headers such as "Results and Discussion" or "Results & Discussion"
       are normalised and mapped to the canonical name "results".
    This approach avoids a heavy multi-line regex and therefore uses less CPU
    and memory on large PDFs.
    """
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

    # Map merged variations to results
    merged_variants = [
        'results and discussion',
        'results & discussion',
        'results discussion'
    ]

    sections: dict[str, str] = {}
    current_header = 'preamble'
    buffer_lines: list[str] = []

    def flush():
        if buffer_lines:
            sections[current_header] = '\n'.join(buffer_lines).strip()
            buffer_lines.clear()

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            buffer_lines.append(raw_line)
            continue
        # Normalise the line for header comparison
        norm = re.sub(r'[^a-z& ]', '', line.lower())
        norm = re.sub(r'\s+', ' ', norm).strip()
        if norm in merged_variants:
            flush()
            current_header = 'results'  # treat merged header as results
            continue
        if norm in known_headers:
            flush()
            current_header = known_headers[norm]
            continue
        # ALL CAPS heuristic (section header without punctuation)
        if line.isupper() and len(line) >= 3 and len(line.split()) <= 6:
            flush()
            current_header = re.sub(r'\s+', '', line.lower())
            continue
        buffer_lines.append(raw_line)
    # flush the last section
    flush()

    # If results missing but discussion exists (merged later), copy discussion
    if 'results' not in sections and 'discussion' in sections:
        sections['results'] = sections['discussion']

    return sections
    # Build a regex pattern to match section headers with arbitrary spaces and case
    section_names = [
        'abstract',
        'results',
        'material and methods',
        'keywords',
        # merged headers
        'results and discussion',
        'results & discussion',
        'results discussion'
    ]
    patterns = [spaced_re(s) for s in section_names]
    # Join all patterns with | for alternation
    header_pattern = '(' + '|'.join(patterns) + ')'
    # Section ends at a line that looks like a new header or end of string
    pattern = header_pattern + r'(.*?)(?=\n[A-Z][A-Za-z ]{3,}:|\Z)'
    sections = {}
    for match in re.finditer(pattern, text, re.IGNORECASE | re.DOTALL):
        # Normalize header by removing spaces and lowering case
        header_raw = re.sub(r'\s+', ' ', match.group(1)).strip().lower()
        header = header_raw.replace('&', 'and').replace('  ', ' ')
        header = re.sub(r'[^a-z ]', '', header)
        # Map merged results+discussion header to 'results'
        if header.startswith('results') and 'discussion' in header:
            header = 'results'
        sections[header] = match.group(2)
    # --- Add discussion section if not present, by normalized word ---
    # Ensure discussion section is detected or inferred
    if 'discussion' not in sections:
        # Normalize text for search
        norm_text = re.sub(r'\s+', '', text).lower()
        idx = norm_text.find('discussion')
        if idx != -1:
            # Find corresponding index in original text
            # Go through original text, skipping whitespace, to find the char index
            count = 0
            orig_idx = None
            for i, c in enumerate(text):
                if not c.isspace():
                    if count == idx:
                        orig_idx = i
                        break
                    count += 1
            if orig_idx is not None:
                sections['discussion'] = text[orig_idx:]
    # If results section is still missing but discussion exists, treat discussion as results (merged section)
    if 'results' not in sections and 'discussion' in sections:
        sections['results'] = sections['discussion']
    return sections

def extract_insect_species(section_text, insect_db, debug_pdf=None):
    found = set()
    candidates = []
    genus_species_in_db = set()  # Genus for which a genus-species pair is present in the DB
    genus_species_matches = set()
    # Find genus-species pairs
    for match in re.findall(r'\b[A-Z][a-z]{2,}\s*[ \n\r\t]+\s*[a-z]{3,}\b', section_text):
        norm_match = normalize_text(match)
        norm_match = normalize_insect_name(norm_match)
        candidates.append((match, norm_match))
        if norm_match in insect_db:
            genus_species_matches.add(match)
            genus_species_in_db.add(norm_match.split()[0])
    # Find abbreviated genus-species pairs, e.g. D. melanogaster
    abbrev_matches = re.findall(r'\b([A-Z])\.\s*([a-z]{3,})\b', section_text)
    # Build genus initial -> genus name map from taxon DB
    genus_map = {}
    for name in insect_db:
        genus = name.split()[0] if ' ' in name else name
        if genus and len(genus) > 0:
            genus_map.setdefault(genus[0], set()).add(genus)
    abbrev_matches_found = set()
    for initial, species in abbrev_matches:
        initial = initial.lower()
        for genus in genus_map.get(initial, []):
            norm_abbrev = f"{genus} {species}"
            if norm_abbrev in insect_db:
                abbrev_matches_found.add(f"{genus} {species}")
                genus_species_in_db.add(genus)
                candidates.append((f"{initial.upper()}. {species}", norm_abbrev))
    # Find single-word genus matches (capitalized, at least 3 letters)
    genus_only_matches = set()
    for genus in re.findall(r'\b[A-Z][a-z]{2,}\b', section_text):
        genus_norm = normalize_text(genus)
        genus_norm = normalize_insect_name(genus_norm)
        # Only add genus if no genus-species pair for this genus is present in the DB
        if genus_norm and genus_norm in insect_db and genus_norm not in genus_species_in_db:
            genus_only_matches.add(genus)
    if debug_pdf:
        logging.info(f"DEBUG: Candidates in {debug_pdf}: {candidates}")
        logging.info(f"DEBUG: insect_db sample: {list(sorted(insect_db))[:50]}")
    # Prefer full genus-species matches
    if genus_species_matches:
        return list(genus_species_matches)
    elif abbrev_matches_found:
        return list(abbrev_matches_found)
    elif genus_only_matches:
        # Only return genus-species pairs that were actually matched in this PDF (from this extraction)
        # This means: only genus-species pairs found in this run, not all from the DB
        # If no genus-species or abbrev found, genus_only_matches are returned as-is
        return []
    else:
        return []


_stress_re = re.compile(r'\b(?:\w+\s)?(?:stress|tolerance)\b', re.IGNORECASE)
# List of common/irrelevant words to exclude from stress sources
COMMON_STRESS_WORDS = set([
    'eukaryotic','certain','s','anol','against','different','time','this','mass','produce','fects','tive','induced','cross','facilitated','factors','better','high','higher','strong','weak','under','vation','increased','decreased','tion','its','it','itself','have','has','had','is','are','was','were','be','been','being','do','does','did','will','shall','can','could','may','might','must','should','the', 'a', 'in', 'to', 'of', 'for', 'by', 'with', 'as', 'on', 'at', 'from', 'and', 'or', 'nor', 'but', 'so', 'yet', 'all', 'any', 'during', 'meets', 'multiple', 'post', 'subsequent', 'when', 'adverse', 'affects', 'all', 'although', 'among', 'analysis', 'avoid', 'broadens', 'causes', 'causing', 'confers', 'controlling', 'critical', 'dative', 'different', 'greater', 'hardening', 'ical', 'idative', 'improved', 'induce', 'induces', 'lution', 'mal', 'mediate', 'otic', 'pathways', 'severe', 'showed', 'stronger', 'sublethal', 'term', 'their', 'that', 'upon', 'various', 'lution', 'avoid'
])

def extract_stress_sources(text):
    # Match 'oxidative stress', 'cold tolerance', etc. (one word before or just the keyword)
    matches = _stress_re.findall(text)
    # Remove 'tolerance' and 'stress', normalize and deduplicate
    sources = set()
    for m in matches:
        cleaned = m.lower().replace('tolerance', '').replace('stress', '')
        cleaned = cleaned.strip().replace('-', ' ').replace('_', ' ')
        cleaned = ' '.join(cleaned.split())  # collapse multiple spaces
        if cleaned:
            sources.add(cleaned)
    # Remove common/irrelevant words
    filtered_sources = [s for s in sources if s not in COMMON_STRESS_WORDS]
    # Deduplicate and return as sorted list
    return sorted(filtered_sources)


from functools import lru_cache

# Cache extracted text for debug PDFs (by path)
@lru_cache(maxsize=128)
def _cached_extract_text(pdf_path):
    return extract_text(pdf_path)

def process_pdf(pdf_path, hmdb_db, insect_db, genus_to_species_counter):
    pdf_name = os.path.basename(pdf_path)
    # Use cache for debug PDFs
    if pdf_name in debug_pdfs:
        text = _cached_extract_text(pdf_path)
    else:
        text = extract_text(pdf_path)
    sections = extract_sections(text)
    # Use the first 10 non-empty lines as the title block
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    title_block = ' '.join(lines[:10])
    # Build relevant text for metabolite extraction.
    metabolite_sections = [
        title_block,
        sections.get('abstract', ''),
        sections.get('results', ''),
        sections.get('materialsandmethods', ''),
        sections.get('keywords', ''),
        sections.get('introduction', ''),
        sections.get('discussion', ''),
        sections.get('conclusion', ''),
        sections.get('references', '')
    ]
    relevant_metabolite_text = ' '.join(metabolite_sections)
    # Build relevant text for insect species extraction (exclude discussion, conclusion)
    insect_sections = [
        title_block,
        sections.get('abstract', ''),
        sections.get('results', ''),
        sections.get('materialsandmethods', ''),
        sections.get('keywords', ''),
        sections.get('introduction', '')
    ]
    relevant_insect_text = ' '.join(insect_sections)
    # Extract tokens and generate n-grams for metabolite search
    tokens = extract_metabolite_tokens(relevant_metabolite_text)
    ngrams = generate_ngrams(tokens, max_n=8)
    # Normalize n-grams and match against HMDB
    norm_ngrams = set(normalize(ng) for ng in ngrams)
    # Log per-PDF diagnostics
    logging.info(f"PDF {pdf_name}: Tokens={len(tokens)}, N-grams={len(ngrams)}, Unique normalized n-grams={len(norm_ngrams)}")
    if pdf_name in debug_pdfs:
        print(f"[DEBUG] Total tokens: {len(tokens)} | Total n-grams: {len(ngrams)} | Unique normalized n-grams: {len(norm_ngrams)}")
        print(f"[DEBUG] Sample tokens: {tokens[:30]}")
        print(f"[DEBUG] Sample n-grams: {ngrams[:30]}")
    
    # Exact matching
    matched = [hmdb_db[ng] for ng in norm_ngrams if ng in hmdb_db]
    matched_norm_ngrams = set(ng for ng in norm_ngrams if ng in hmdb_db)
    
    # Log exact match count
    logging.info(f"PDF {pdf_name}: Exact matched metabolites (before fuzzy)={len(matched)}")
    if pdf_name in debug_pdfs:
        print(f"[DEBUG] Exact matched count: {len(matched)}")
    
    # Fuzzy fallback: for non-matched n-grams, use RapidFuzz token_set_ratio >= 90
    non_matched_ngrams = [ng for ng in norm_ngrams if ng not in matched_norm_ngrams]
    fuzzy_matched = []
    fuzzy_match_count = 0
    
    if non_matched_ngrams:
        # Limit fuzzy search to control false positives - only check longer n-grams
        fuzzy_candidates = [ng for ng in non_matched_ngrams if len(ng.split()) >= 2]
        
        for candidate in fuzzy_candidates:
            best_match = None
            best_score = 0
            # Check against all HMDB keys
            for hmdb_key in hmdb_db.keys():
                score = fuzz.token_set_ratio(candidate, hmdb_key)
                if score >= 90 and score > best_score:
                    best_score = score
                    best_match = hmdb_key
            
            if best_match:
                fuzzy_matched.append(hmdb_db[best_match])
                fuzzy_match_count += 1
    
    # Combine exact and fuzzy matches
    matched.extend(fuzzy_matched)
    
    # Log fuzzy match count
    if fuzzy_match_count > 0:
        logging.info(f"PDF {pdf_name}: Fuzzy matched metabolites={fuzzy_match_count}")
    
    # Log match count
    logging.info(f"PDF {pdf_name}: Total matched metabolites (before single-element filter)={len(matched)}")
    if pdf_name in debug_pdfs:
        print(f"[DEBUG] Total matched count (exact + fuzzy, before single-element filter): {len(matched)}")
    
    matched = [m for m in matched if not is_single_element(m['CHEMICAL_FORMULA'])]
    
    # Log unique metabolites by INCHIKEY or NAME
    unique_metabolites = set()
    for m in matched:
        unique_key = m.get('INCHIKEY') if m.get('INCHIKEY') else m.get('NAME', '')
        if unique_key:
            unique_metabolites.add(unique_key)
    logging.info(f"PDF {pdf_name}: Unique matched metabolites (after filtering)={len(unique_metabolites)}")
    
    if not matched:
        logging.info(f"No metabolites found in {os.path.basename(pdf_path)}; skipping.")
        return None
    # Insect species from allowed sections
    insect_species = extract_insect_species(relevant_insect_text, insect_db, debug_pdf=pdf_name if pdf_name in debug_pdfs else None)
    # If only genus names are found, use the most common species for that genus from the cache
    normalized_insect_species = []
    norm_species = [normalize_insect_name(s) for s in insect_species]
    # If all are single-word (genus only), use genus_to_species_counter
    if all(' ' not in s for s in norm_species):
        for genus in norm_species:
            if genus in genus_to_species_counter and genus_to_species_counter[genus]:
                most_common = genus_to_species_counter[genus].most_common(1)
                normalized_insect_species.extend([f"{genus} {species}" for species, _ in most_common])
    else:
        normalized_insect_species = sorted(set(norm_species))
    if pdf_name in debug_pdfs:
        print(f"[DEBUG] Insect species found: {insect_species}")
        print("[DEBUG] Normalized insect species candidates and DB presence:")
        for s in insect_species:
            norm = normalize_insect_name(s)
            print(f"    '{s}' (normalized: '{norm}') in DB: {norm in insect_db}")
        print(f"[DEBUG] Sample normalized insect DB: {list(insect_db)[:50]}")
    if pdf_name in debug_pdfs:
        print(f"[DEBUG] Insect species found: {insect_species}")
        print("[DEBUG] Normalized insect species candidates and DB presence:")
        for s in insect_species:
            norm = normalize_insect_name(s)
            print(f"    '{s}' (normalized: '{norm}') in DB: {norm in insect_db}")
        print(f"[DEBUG] Sample normalized insect DB: {list(insect_db)[:50]}")
    if not normalized_insect_species:
        logging.info(f"PDF {pdf_name}: Skipped due to missing insect species (gate check failed)")
        print(f"[DEBUG] No insect species found, skipping.")
        return None
    # Stress sources: can be from all sections
    relevant_stress_text = ' '.join([
        title_block,
        sections.get('abstract', ''),
        sections.get('results', ''),
        sections.get('materialsandmethods', ''),
        sections.get('keywords', ''),
        sections.get('introduction', ''),
        sections.get('discussion', ''),
        sections.get('conclusion', ''),
        sections.get('references', '')
    ])
    stress_sources = extract_stress_sources(relevant_stress_text)
    # DOI
    dois = list(dict.fromkeys(extract_doi(text)))  # Preserve order, deduplicate
    doi_reference = dois[0] if dois else ''
    # Compile results
    results = []
    for m in matched:
        logp = get_logp(m['INCHIKEY'])
        superclass = get_npclassifier_superclass(m['SMILES'], m['INCHIKEY'])
        results.append({
            'Metabolite name': m.get('NAME', ''),
            'CHEMICAL_FORMULA': m.get('CHEMICAL_FORMULA', ''),
            'Average Mass': m.get('AVERAGE_MASS', ''),
            'Mono Mass': m.get('MONO_MASS', ''),
            'LogP': logp,
            'superclass': superclass,
            'stress sources': ', '.join(sorted(set(s for s in stress_sources if s))),
            'insect species': ', '.join(normalized_insect_species),
            'DOI references': doi_reference
        })
        if pdf_name in debug_pdfs:
            print(f"[DEBUG] Results: {results}")
    return results

from collections import defaultdict, Counter

def main():
    from collections import defaultdict, Counter
    hmdb_db = load_hmdb(HMDB_PATH)
    insect_db = load_taxon(TAXON_PATH)
    pdf_files = [os.path.join(PDF_DIR, f) for f in os.listdir(PDF_DIR) if f.endswith('.pdf')]
    # First pass: build genus_to_species_counter from all PDFs
    genus_to_species_counter = defaultdict(Counter)
    pdf_to_species = {}
    for pdf in pdf_files:
        text = extract_text(pdf)
        sections = extract_sections(text)
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        title_block = ' '.join(lines[:10])
        insect_sections = [
            title_block,
            sections.get('abstract', ''),
            sections.get('results', ''),
            sections.get('materialsandmethods', ''),
            sections.get('keywords', ''),
            sections.get('introduction', '')
        ]
        relevant_insect_text = ' '.join(insect_sections)
        # Extract all genus-species found in this PDF
        species_found = []
        for match in re.findall(r'\b[A-Z][a-z]{2,}\s*[ \n\r\t]+\s*[a-z]{3,}\b', relevant_insect_text):
            norm_match = normalize_text(match)
            norm_match = normalize_insect_name(norm_match)
            if norm_match in insect_db and ' ' in norm_match:
                genus, species = norm_match.split()[:2]
                genus_to_species_counter[genus][species] += 1
                species_found.append(norm_match)
        # Abbreviated genus-species (e.g. D. melanogaster)
        abbrev_matches = re.findall(r'\b([A-Z])\.\s*([a-z]{3,})\b', relevant_insect_text)
        genus_map = {}
        for name in insect_db:
            genus = name.split()[0] if ' ' in name else name
            if genus and len(genus) > 0:
                genus_map.setdefault(genus[0], set()).add(genus)
        for initial, species in abbrev_matches:
            initial = initial.lower()
            for genus in genus_map.get(initial, []):
                norm_abbrev = f"{genus} {species}"
                if norm_abbrev in insect_db:
                    genus_to_species_counter[genus][species] += 1
                    species_found.append(norm_abbrev)
        pdf_to_species[os.path.basename(pdf)] = set(species_found)
    # Second pass: process each PDF with full genus_to_species_counter
    all_results = []
    for pdf in pdf_files:
        try:
            res = process_pdf(pdf, hmdb_db, insect_db, genus_to_species_counter)
            if res:
                all_results.extend(res)
        except Exception as e:
            logging.exception(f"Error processing {pdf}: {e}")
    if all_results:
        # Merge duplicates by InChIKey (or Metabolite name if InChIKey missing)
        merged = {}
        for row in all_results:
            key = row.get('INCHIKEY') or row.get('Metabolite name')
            if key not in merged:
                merged[key] = row.copy()
            else:
                # Merge unique values for each field
                for k, v in row.items():
                    if k in ['stress sources', 'insect species', 'DOI references']:
                        old = set(merged[key][k].split(', ')) if merged[key][k] else set()
                        new = set(v.split(', ')) if v else set()
                        # For stress sources, remove empty, normalize and deduplicate
                        merged[key][k] = ', '.join(sorted(s.strip() for s in (old | new) if s.strip()))
        df = pd.DataFrame(merged.values())
        df.to_excel(OUTPUT_XLSX, index=False)
        print(f"Results saved to {OUTPUT_XLSX}")
    else:
        print("No metabolites extracted from any PDFs.")

if __name__ == '__main__':
    main()


