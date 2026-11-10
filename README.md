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
normalization); the displayed query terms retain the extracted spelling.
`Matched name/synonym` shows the exact name or synonym from the source database
that matched, retaining its original spelling. For an `EXACT` match, it lists
all matched database terms, one per line. For a `FUZZY` match, `Fuzzy query`
contains the normalized query and `Fuzzy match` contains the normalized
database term. The report omits HMDB/LMP IDs; use `Match type`, `Database`, and
`Matched name/synonym` to review the match. A database match does not by itself
establish that the identification is valid.
