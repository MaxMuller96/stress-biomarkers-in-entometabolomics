# stress-biomarkers-in-entometabolomics
An algorithm that uses natural language processing for listing of stress biomarkers in insects in metabolomics derived from published papers.
The algorithm is written in Python in file stressbiomarkers_insects.py.

The referenced papers need to be downloaded manually.
Databases need to be downloaded locally from GBIF.org for the taxa and HMDB.ca for the metabolites.

**Startup**
Install python 3.14 or lower and an IDE like Spyder or Jupyter.
Download requirements.txt and place it in the project directory.
In the terminal execute: cd "directory" and py -m pip install -r requirements.txt
If the packages are installed, open the code in Spyder or Jupyter.
Run the code.

In the `metabolite_hits` sheet, `Fuzzy query` records the text used to
identify each metabolite. For an `EXACT` match it lists all matching names,
synonyms, systematic names, or abbreviations found in the paper, one per line.
Matching is exact after normalization (including case and Greek-letter
normalization); the displayed terms retain the extracted spelling. For a
`FUZZY` match, the column contains the normalized query and `Fuzzy match`
contains the normalized database term. Use `Match type`, `Database`, and
`HMDB ID/LMP` to review the match; a database match does not by itself establish
that the identification is valid.
