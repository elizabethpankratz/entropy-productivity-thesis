# -*- coding: utf-8 -*-
# Extracts candidate derivations matching the queries in to_sample.csv up to sample sizes of 100,000.
# Saves each as its own file in raw_samples/.

from SeaCOW import Query, Nonprocessor, ConcordanceLoader
import pandas as pd

CORPUS      = 'decow16b'

morph_df    = pd.read_csv('queries.csv')
MORPHS      = morph_df['morph'].dropna()
QUERIES     = morph_df['query'].dropna()
NUM_MORPHS  = len(MORPHS)


# ======================================================


def conduct_query(cql_string, corpus):
  """
  Queries the given corpus.
  
  Args:
    cql_string: string in CQL format
    corpus: string representing the corpus to query ('decow16a-nano', 'decow16b', 'encow16a', 'encow16a-nano')
    subcorpus: string representing the subcorpus to query 
  Returns:
    List of dicts containing the concordances, of length max_hits
  """
  
  q = Query()
  q.corpus          = corpus  
  q.string          = cql_string
  q.max_hits        = 100000
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

# Draw samples  from DECOW for each suffix.
for m_idx in range(NUM_MORPHS):
  curr_morph = MORPHS[m_idx]
  curr_query = QUERIES[m_idx]

  concs = conduct_query(curr_query, CORPUS)

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

  # Wherever 'compana' isn't '_' (i.e., wherever it's a compound), split on _ and extract first and last elements.
  conc_df['cpd.N1'] = pd.np.where(conc_df.compana != '_',
                                conc_df.compana.str.split("_").str[0],
                                '')

  conc_df['cpd.N2'] = pd.np.where(conc_df.compana != '_',
                                conc_df.compana.str.split("_").str[-1],
                                '')
                                
  # Add morpheme as column to df.
  conc_df['morph'] = [curr_morph] * len(conc_df)

  # Rearrange columns.
  conc_df = conc_df[['morph', 'doc.url','doc.id','s.idx','word','lemma', 'compana', 'cpd.N1', 'cpd.N2']]

  # Save as CSV.
  conc_df.to_csv('1_raw_samples/%s.csv' % curr_morph, index=False, encoding='UTF-8')

  print 'Done %s/%s: %s' % (m_idx+1, NUM_MORPHS, curr_morph)

