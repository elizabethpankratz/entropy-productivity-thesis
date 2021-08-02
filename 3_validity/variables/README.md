
**Script:**
- `compute_variables.ipynb`: Computes the three factors to be used in the model: frequency ratio, semantic relatedness, and junctural phonotactics.
  - In:
    - Frequency ratio: contents of `../../1_data/35_samples/6_backform_base_cutoff/`
    - Semantic relatedness: `backformer_two.py`, `DErivBase-v2.0-probabilities.txt`
    - Junctural phonotactics: `../simplexes/junc_data/junctures_tokenbased.csv`,  contents of `../../1_data/35_samples/7_analysis_samples/`
    - Entropy: contents of `../../1_data/35_samples/7_analysis_samples/`
  - Out:
    - `sfx_data.csv`

**Module:**
- `backformer_two.py`: Version 2 of `backformer` module, now updated based on rules that were discovered to be missing while annotating the generated bases.

**Data files:**
- `DErivBase-v2.0-probabilities.txt`: From DErivBase 2.0, the learned probabilities that each pair of words is semantically related.
- `sfx_data.csv`: The dataframe at the heart of the analysis; all productivity factors and entropy for each suffix.

