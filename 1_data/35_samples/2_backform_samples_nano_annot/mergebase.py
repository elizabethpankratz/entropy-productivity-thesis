# Merges the annotations in the given X_unique_annot.csv file with the lemmas in X_base_annot.csv.
# Requires command line argument of suffix to process, e.g.:
# > python mergebase.py -ung


import pandas as pd
import sys

SFX = sys.argv[1]

df_uni = pd.read_csv(SFX + '_unique_annot.csv')
df_all = pd.read_csv(SFX + '_base_annot.csv')

df_merged = df_all.drop(columns=['true_lemma', 'merge']).merge(df_uni, on='lemma')
df_merged = df_merged[['lemma', 'unique_candidates', 'pos', 'lemma_freq', 'base_freq', 'true_lemma', 'true_base', 'query_by_hand', 'query_pos', 'merge']]

df_merged.to_csv(SFX + '_base_annot.csv', index=False)

print('Merged lemma annotations for', SFX)
