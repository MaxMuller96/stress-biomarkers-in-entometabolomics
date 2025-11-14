# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 14:24:33 2025

@author: maxmu
"""

import os
import re
import csv
import logging
import pandas as pd
from tqdm import tqdm
from pdfminer.high_level import extract_text
import requests

# Paths (update as needed)
PDF_DIR = r".../Stressbiomarkers_insects/pdfs"
HMDB_PATH = r".../hmdb_metabolites/metabolites-2025-06-20.txt"
TAXON_PATH = r".../Stressbiomarkers_insects/Taxon/Taxon.txt"
LOG_PATH = "skipped_papers.log"
OUTPUT_XLSX = "metabolite_results.xlsx"

# Setup logging
logging.basicConfig(filename=LOG_PATH, level=logging.INFO, format='%(asctime)s %(message)s')

def normalize(text):
    # Remove hyphens, delimiters, spaces, dashes, etc., and lowercase
    return re.sub(r'[-–—,;:\s]', '', text.lower())

def extract_words(text):
    # Remove hyphens at line breaks, normalize, and extract words
    text = re.sub(r'-\s*\n\s*', '', text)
    words = re.findall(r'\b\w+\b', text)
    return [w for w in words if len(w) >= 3]

def load_hmdb(path):
    db = {}
    with open(path, encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
            name_norm = normalize(row['name'])
            db[name_norm] = row
    return db

def load_taxon(path):
    insects = set()
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line:
                insects.add(line.lower())
    return insects

def get_logp(inchikey):
    # PubChem REST API
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/XLogP/JSON'
    try:
        r = requests.get(url, timeout=10)
        if r.ok:
            props = r.json()['PropertyTable']['Properties'][0]
            return props.get('XLogP', None)
    except Exception:
        pass
    return None

def get_npclassifier_superclass(smiles, inchikey):
    # NPClassifier API (example endpoint)
    url = 'https://npclassifier.ucsd.edu/classify'
    try:
        r = requests.post(url, json={'smiles': smiles or '', 'inchikey': inchikey or ''}, timeout=10)
        if r.ok:
            data = r.json()
            return data.get('superclass', None)
    except Exception:
        pass
    return None

def extract_sections(text):
    # Extract sections: results, abstract, materials and methods
    pattern = r'(abstract|results|materials and methods)(.*?)(?=\n[A-Z][A-Za-z ]{3,}:|\Z)'  # crude section splitter
    sections = {}
    for match in re.finditer(pattern, text, re.IGNORECASE | re.DOTALL):
        sections[match.group(1).lower()] = match.group(2)
    return sections

def extract_insect_species(section_text, insect_db):
    # Regex: \b[A-Z][a-z]{2,} [a-z]{3,}\b
    found = set()
    for match in re.findall(r'\b[A-Z][a-z]{2,} [a-z]{3,}\b', section_text):
        if match.lower() in insect_db:
            found.add(match)
    return list(found)

def extract_stress_sources(text):
    # Regex: \b(?:\w+ )?stress\b
    return list(set(re.findall(r'\b(?:\w+ )?stress\b', text, re.IGNORECASE)))

def extract_doi(text):
    # Regex: \b10\.\d{4,9}/[-._;()/:A-Z0-9]+\b
    return re.findall(r'\b10\.\d{4,9}/[-._;()/:A-Z0-9]+\b', text, re.IGNORECASE)

def is_single_element(formula):
    # Only 1 element, e.g. 'Cl', 'C', 'O2' (<=2 letters)
    return len(re.findall(r'[A-Z][a-z]?', formula)) == 1

def process_pdf(pdf_path, hmdb_db, insect_db):
    text = extract_text(pdf_path)
    words = extract_words(text)
    norm_words = set(normalize(w) for w in words)
    # Match metabolites
    matched = [hmdb_db[w] for w in norm_words if w in hmdb_db]
    matched = [m for m in matched if not is_single_element(m['Chemical_Formula'])]
    if not matched:
        logging.info(f"No metabolites found in {os.path.basename(pdf_path)}")
        return None
    # Extract sections
    sections = extract_sections(text)
    relevant_text = ' '.join(sections.get(s, '') for s in ['abstract', 'results', 'materials and methods'])
    # Insect species
    insect_species = extract_insect_species(relevant_text, insect_db)
    if not insect_species:
        logging.info(f"No insect species found in {os.path.basename(pdf_path)}; metabolites removed")
        return None
    # Stress sources
    stress_sources = extract_stress_sources(relevant_text)
    # DOI
    dois = extract_doi(text)
    # Compile results
    results = []
    for m in matched:
        logp = get_logp(m['INCHIKEY'])
        superclass = get_npclassifier_superclass(m['SMILES'], m['INCHIKEY'])
        results.append({
            'Metabolite name': m['name'],
            'Chemical formula': m['Chemical_Formula'],
            'Average Mass': m['Average_Mass'],
            'Mono Mass': m['Mono_Mass'],
            'LogP': logp,
            'superclass': superclass,
            'stress sources': ', '.join(stress_sources),
            'insect species': ', '.join(insect_species),
            'DOI references': ', '.join(dois),
        })
    return results

def main():
    hmdb_db = load_hmdb(HMDB_PATH)
    insect_db = load_taxon(TAXON_PATH)
    all_results = []
    pdf_files = [os.path.join(PDF_DIR, f) for f in os.listdir(PDF_DIR) if f.lower().endswith('.pdf')]
    for pdf in tqdm(pdf_files, desc="Processing PDFs"):
        try:
            res = process_pdf(pdf, hmdb_db, insect_db)
            if res:
                all_results.extend(res)
        except Exception as e:
            logging.exception(f"Error processing {pdf}: {e}")
    if all_results:
        df = pd.DataFrame(all_results)
        df.to_excel(OUTPUT_XLSX, index=False)
        print(f"Results saved to {OUTPUT_XLSX}")
    else:
        print("No metabolites extracted from any PDFs.")

if __name__ == '__main__':
    main()

