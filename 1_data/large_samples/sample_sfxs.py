# -*- coding: utf-8 -*-
# Extracts all derivations with the given suffixes from DECOW16A-NANO.

from SeaCOW import Query, Nonprocessor, ConcordanceLoader
import pandas as pd

CORPUS        = 'decow16a-nano'  
MATRIX_Q      = '[lemma=".*%s$" & tag="NN"] within <s/>'
SFXS           = ['(h|k)eit', 'schaft', 'nis']
SFXS_PLAIN     = ['heit', 'schaft', 'nis']

# ======================================================


def conduct_query(cql_string, corpus):
  """
  Queries the given corpus.
  
  Args:
    cql_string: string in CQL format
    corpus: string representing the corpus to query ('decow16a-nano', 'decow16b', 'encow16a', 'encow16a-nano')
  Returns:
    List of dicts containing the concordances, of length max_hits
  """
  
  q = Query()
  q.corpus          = corpus  
  q.string          = cql_string
  q.max_hits        = -1
  q.attributes      = ['word', 'lemma', 'compana']
  q.structures      = ['s']
  q.references      = ['doc.url', 'doc.id', 's.idx']
  q.container       = 's'
  q.set_deduplication()
  
  p                 = ConcordanceLoader()
  p.full_structure  = True     # Convert token attributes to dicts as well, otherwise |-separated. 
  q.processor       = p
  q.run()
  
  return p.concordance
  

# ======================================================

# Draw samples for each suffix from DECOW.
for sfx_idx in range(len(SFXS)):
  sfx = SFXS[sfx_idx]
  sfx_plain = SFXS_PLAIN[sfx_idx]

  concs = conduct_query(MATRIX_Q % sfx, CORPUS)

  # Extract the information from each match that we probably care about.
  conc_data = [
    {
      'word': conc['match'][0]['word'],
      'lemma': conc['match'][0]['lemma'],
      'compana': conc['match'][0]['compana'],
      'doc.id': conc['meta']['doc.id'],
      'doc.url': conc['meta']['doc.url'],
      's.idx': conc['meta']['s.idx']
    } for conc in concs]

  # Convert list of dicts to df.
  conc_df = pd.DataFrame(conc_data)
  
  # Wherever 'compana' isn't '_' (i.e., wherever it's a compound), split on _ and extract final element (the head).
  # Wherever it isn't a compound, just keep the value in lemma. This is the new value of the lemma column. 
  # Drop compana column.
  conc_df['lemma'] = pd.np.where(conc_df.compana != '_', 
                                conc_df.compana.str.split("_").str[-1], 
                                conc_df.lemma)
  conc_df.drop(['compana'], axis=1, inplace=True)
  
  # If lemma contains - (e.g. 'DSL-Geschwindigkeit'), split at - and take final element (should
  # always be the derivation).
  conc_df['lemma'] = pd.np.where(conc_df.lemma.str.contains('-'), 
                                conc_df.lemma.str.split("-").str[-1], 
                                conc_df.lemma)
  # Rearrange columns.
  conc_df = conc_df[['doc.url','doc.id','s.idx','word','lemma']]
  
  # Save as CSV.
  conc_df.to_csv('%s.csv' % sfx_plain, index=False, encoding='UTF-8')
  
  print 'Done', sfx
  
