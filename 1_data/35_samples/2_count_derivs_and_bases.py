# -*- coding: utf-8 -*-
# Count the number of times each derivation and each reconstructed base appears in the corpus.
# To be run after 1_sample_sfxs.py.

from SeaCOW import Query, Nonprocessor
import backformer_one as b
import pandas as pd
import numpy as np

CORPUS        = 'decow16a-nano'
SAMPLE_SIZE   = 20000
SFXS          = [
  '-age', '-and', '-ant', '-anz', '-ation',
  '-atur', '-ement', '-end', '-ent', '-enz',
  '-eur', '-iment',
  '-iteur', '-itur', '-ament',
  '-ateur', '-ator',
  '-el',
  '-er',
  '-heit', '-ie', '-ik', '-iker', '-ikum',
  '-ismus', '-ist',
  '-itaet',
  '-ition', '-ium',
  '-ling', '-nis', '-schaft',
  '-ung',
  '-e' # gives UnicodeDecodeError but still saves
  ]


# ======================================================


def get_count(cql, corpus):
  """
  Counts the number of tokens for the given query in the given corpus. Returns number of hits (i.e., tokens).
  """
  
  # Create a Query object and set whatever needs to be set.
  q = Query()
  q.corpus          = corpus                # Lower-case name of the corpusto use.
  q.string          = cql                     # A normal CQL string as used in NoSketchEngine.
  q.max_hits        = -1
  q.attributes      = []                      # For counting, you don't need word attributes.
  q.structures      = []                      # ... you don't need structural attributes.
  q.references      = []                      # ... you don't need reference attrs.
  q.container       = 's'                       # Which container structure should be used? None is OK
                                              # only if class is Nonprocessor.
  
  # Using the deduplicator would NOT change the outcome. Switch off.
  q.set_deduplication(off = True)
  
  # Create a Processor object and attach it to the Query object.
  # The Nonprocessor processor does nothing. You can work with the results
  # yourself in the finalise method or just get the hits value from the
  # query object. It is the concordance as seported by Manatee.
  p                 = Nonprocessor()  # Create a processor object of appropriate type.
  q.processor       = p               # Attach the processor to the query.
  q.run()                             # Run the query.
  return q.hits


# ======================================================

for sfx in SFXS:
  
  # Read in the sample for the current suffix (query was done in 1_sample_sfxs.py).
  curr_sample = pd.read_csv('raw_samples/'+sfx+'.csv')#, encoding='UTF-8')
  
  # Take a random sample of size SAMPLE_SIZE. This will be the basis of the analysed data.
  # (It's way too slow to backform and test bases for 100k tokens, and we know that entropy
  # should stabilise upward of ~10k tokens.)
  # Crashes if the size of curr_sample is < SAMPLE_SIZE, so if that happens, just take the whole sample as sample_subset.
  try:
    rd_idcs = np.random.choice(curr_sample.index, size = SAMPLE_SIZE, replace = False)
    sample_subset = curr_sample.iloc[rd_idcs].reset_index(drop=True)
  except:
    sample_subset = curr_sample
    
  # Save the current sample in random_subsamples/. This sample will form the basis of the rest of the analyses.
  sample_subset.to_csv('2_random_subsamples/' + sfx + '_subsample.csv', index=False, encoding='UTF-8')
  
  # Get the CQL queries for candidate bases for the current suffix via Backformer.
  curr_bases = b.get_bases(sample_subset, sfx)
  curr_cql_df = b.get_cql_from_bases(curr_bases)

  # Conduct queries and add as new column to curr_cql_df.
  curr_cql_df['base_freq'] = [get_count(query, CORPUS) for query in curr_cql_df['cql']]

  # Get frequencies of derivations too (uniquify them and then merge freqs with curr_cql_df).
  lemma_freqs = [{'lemma':lemma, 'lemma_freq':get_count('[lemma="%s" & tag="NN"] within <s/>' % lemma, CORPUS)} for lemma in curr_cql_df['lemma'].unique()]
  curr_cql_df = curr_cql_df.merge( pd.DataFrame(lemma_freqs), on='lemma' )

  # Reorder columns (and drop CQL, don't need it anymore) and save to CSV.
  curr_cql_df = curr_cql_df[['lemma', 'unique_candidates', 'pos', 'lemma_freq', 'base_freq']]
  curr_cql_df.to_csv('2_backform_samples_nano/' + sfx + '_base_freqs.csv', index=False, encoding='UTF-8')

  print 'Done', sfx

