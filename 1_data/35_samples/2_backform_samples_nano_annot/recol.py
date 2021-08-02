import pandas as pd
import os

files = [fn for fn in os.listdir() if fn[-14:] == 'base_annot.csv']
# print(files)

for f in files[:6]:
	df = pd.read_csv(f)
	origcols = df.columns.tolist()
	
	print(origcols)
    #df = df.reindex(origcols[:1] + ['manual_lemma'] + origcols[1:], axis=1)
    #df.to_csv(f, index=False)
    # df[['lemma', 'unique_candidates', 'merge', 'pos', 'lemma_freq', 'base_freq', 'true_lemma', 'true_base', 'query_by_hand']].to_csv(f[:-5] + '.csv', index=False)
