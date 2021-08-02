# Gets unique bases from given _annot.csv file, saves as new document for more efficient annotation.
# This new document will be merged with the _annot.csv file in mergebase.py.
# Requires command line argument of suffix to process, e.g.:
# > python unibase.py -ung

import pandas as pd 
import sys

SFX = sys.argv[1]

FN_ALL = SFX + '_base_annot.csv'
FN_UNI = SFX + '_unique_annot.csv'

df = pd.read_csv(FN_ALL)
df[['lemma', 'true_lemma', 'merge']].drop_duplicates().to_csv(FN_UNI, index=False)

print('Created file', FN_UNI)
