The contents of this directory generate a list of highly common German simplexes.
These will be used for computing the probability of the transitions at morpheme boundaries in the derived words.

# Contents

**Directories:**
- `infiles/`: Files read into the scripts below.
- `outfiles/`: Files produced by the scripts below.

**Scripts:**
- `1_id_probable_simplexes.ipynb`: Reads in most frequent lemmas from DECOW, selects those most likely to be monomorphemic. Preps them for input to SMOR (which happens offstage in between this step and the next).
  - in: `infiles/decow16bx.lp`
  - out: `outfiles/simplex_filtered1.csv`, `outfiles/simplex_filtered1.tosmor`
- `2_filter_w_smor.ipynb`: Reads in SMOR analyses and use them to weed out unwanted lemmas (compounds, numerals, abbrevations, etc.). Saves list of simplexes for manual annotation.
  - in: `infiles/simplex_filtered1.smored`, `infiles/simplex_filtered1.csv`
  - out: `outfiles/simplex_filtered2.csv`, `outfiles/simplex_filtered2_toannot.csv`
- `3_merge_annots.ipynb`: Merges manual annotation with the pruned data, yielding a file that contains mainly simplexes (and probably some unclear borderline cases), their POS, and their frequency in DECOW16B.
  - in: `outfiles/simplex_filtered2.csv`, `infiles/simples_filtered2_annotated.csv`
  - out: `outfiles/simplex_filtered3.csv`
- `4_count_juncture_freq.ipynb`: Computes the frequency of all bigraphs in the simplexes found in Steps 1--3.
  - in: `outfiles/simplex_filtered3.csv`
  - out: `junc_tokenbased.csv`.

**Data files:**
- `junctures_tokenbased.csv`: The token-based probability that each bigraph appears in German simplexes.
