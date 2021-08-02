This material appears in Section 2 as well as in Appendix A.

**Directories:**
- `iterdata/`: All bootstrapped samples for *-heit*, *-nis*, and *-schaft*, with and without replacement, saved as frequency distributions; also, misc. files generated in `bootstrap_prod_measures.Rmd` that would be too costly to generate every time.
- `imgs/`: All figures used in Section 2.
- `bootstrap_prod_measures_files`, `other_plots_files`: Supporting files for markdown.

**Scripts:**
- `gen_bootstrap_samples.ipynb`
  - In: contents of `../1_data/large_samples/`
  - Out: `iterdata/freqdist_iter_500.csv` and `iterdata/freqdist_iter_500_wrepl.csv`
- `bootstrap_prod_measures.Rmd`
  - In: contents of `iterdata/`
  - Out: most contents of `iterdata/`, most contents of `imgs/`
- `other_plots.Rmd`
  - In: nothing
  - Out: `imgs/freqdist.pdf` and `imgs/zipf-selfsimil.pdf`
