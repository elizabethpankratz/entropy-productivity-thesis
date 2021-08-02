# -*- coding: utf-8 -*-

import pandas as pd

def prep_sfx_df(df_raw):
    """
    If there are compounds in the data, then replace the value of 'lemma' with the compound head, in 'cpd.N2'.
    If lemma contains -, split at - and take final element (also the head).
    Drop any rows where the lemma contains |; this is a case of ambiguous lemmatisation.
    
    Arg:
        df_raw: pandas df containing lemma, compana, cpd.N1, cpd.N2 columns, straight out of SeaCOW
    Returns:
        df with cpd.N2 integrated into lemma and cols compana and cpd.N1 dropped.
    """
    df_in = df_raw.copy()
    df_in['lemma'] = pd.np.where(df_raw['cpd.N2'].isnull(),
                      df_raw['lemma'],
                      df_raw['cpd.N2'])
    df_in['lemma'] = pd.np.where(df_in.lemma.str.contains('-'),
                                df_in.lemma.str.split('-').str[-1],
                                df_in.lemma)
    df_in.drop(['compana', 'cpd.N1', 'cpd.N2'], axis=1, inplace=True)
    df_in = df_in[~df_in['lemma'].str.contains('|', regex=False)].reset_index(drop=True)
    
    return df_in
    
def prep_pfx_df(df_raw):
    """
    If there are compounds in the data, then replace the value of 'lemma' with the compound N1, in 'cpd.N1'.
    Drop any rows where the lemma contains |; this is a case of ambiguous lemmatisation.

    Arg:
        df_raw: pandas df containing lemma, compana, cpd.N1, cpd.N2 columns, straight out of SeaCOW
    Returns:
        df with cpd.N1 integrated into lemma and cols compana and cpd.N1 dropped.
    """
    df_in = df_raw.copy()
    df_in['lemma'] = pd.np.where(df_raw['cpd.N1'].isna(),
                      df_raw['lemma'],
                      df_raw['cpd.N1'])
    df_in.drop(columns=['compana', 'cpd.N1', 'cpd.N2'], inplace=True)
    df_in = df_in[~df_in['lemma'].str.contains('|', regex=False)].reset_index(drop=True)

    return df_in

def prep_circfx_df(df_raw):
    """
    If there are compounds in the data, then drop them. Not an example of a derivation with Ge-xxx-e or Ge-xxx-sel.
    Drop any rows where the lemma contains |; this is a case of ambiguous lemmatisation.
    
    Arg:
        df_raw: pandas df containing lemma, compana, cpd.N1, cpd.N2 columns, straight out of SeaCOW
    Returns:
        df with rows containing a compound analysis dropped, and cols compana, cpd.N1, cpd.N2 dropped.
    """
    df_in = df_raw.copy()
    df_in = df_in[df_in['cpd.N1'].isna()].drop(columns=['compana', 'cpd.N1', 'cpd.N2'])
    df_in = df_in[~df_in['lemma'].str.contains('|', regex=False)].reset_index(drop=True)
    
    return df_in

def get_ament_ateur_bases(df_in):
    """
    Backforms the contents of column 'lemma' to the possible stems and saves candidates as new columns.
    (Corresponds to DErivBase dVN24, dVN40)

    Arg:
        df_in: dataframe containing derivation tokens in column 'lemma'
    Returns:
        df with new columns containing potential bases
    """
    df = df_in.copy()
    
    base_v_ieren = []
    base_v_izieren = []
    
    for lemma in df['lemma']:
        base_strip = lemma[:-5].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')
        
        base_ieren = base_strip + 'ieren'
        base_izieren = base_strip[:-1] + 'zieren' if base_strip[-1] == 'k' else ''
        
        base_v_ieren.append(base_ieren)
        base_v_izieren.append(base_izieren)

    df['base_cand1_VVINF'] = base_v_ieren
    df['base_cand2_VVINF'] = base_v_izieren
    
    return df

def get_ator_bases(df_in):
    """
    Backforms the contents of column 'lemma' to the possible stems and saves candidates as new columns.
    (Corresponds to DErivBase dVN20)

    Arg:
        df_in: dataframe containing derivation tokens in column 'lemma'
    Returns:
        df with new columns containing potential bases
    """
    df = df_in.copy()
    
    base_v_ieren = []
    base_v_izieren = []
    
    for lemma in df['lemma']:
        base_strip = lemma[:-4].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')
        
        if len(base_strip) < 1:  # Sometimes we get lemma = "Tor", ignore
            base_ieren = ''
            base_izieren = ''
        else:
            base_ieren = base_strip + 'ieren'
            base_izieren = base_strip[:-1] + 'zieren' if base_strip[-1] == 'k' else ''
        
        base_v_ieren.append(base_ieren)
        base_v_izieren.append(base_izieren)

    df['base_cand1_VVINF'] = base_v_ieren
    df['base_cand2_VVINF'] = base_v_izieren
    
    return df

def get_e_bases(df_in):
    """
    Backforms the contents of column 'lemma' to the possible stems and saves candidates as new columns.
    (Corresponds to DErivBase dAN01, dVN04, dVN08)

    Arg:
        df_in: dataframe containing derivation tokens in column 'lemma'
    Returns:
        df with new columns containing potential bases
    """
    df = df_in.copy()
    
    base_adj = []
    base_adj_nouml = []
    base_v_en = []
    
    irreg = {'gab':'geb', 'sprach':'sprech'}
    
    for lemma in df['lemma']:
        base_strip = lemma[:-1].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')
        base_nouml = base_strip.replace('ü', 'u').replace('ä', 'a').replace('ö', 'o')
        
        # Check if the base is in the dict of irregular stems; if it is, replace it with the value string.
        for key, val in irreg.items():
            if key in base_strip:
                base_strip = base_strip.replace(key, val)
        
        base_adj.append(base_strip.capitalize())
        base_adj_nouml.append(base_nouml if any([uml in base_strip for uml in ['ä', 'ö', 'ü']]) else '')
        base_v_en.append(base_strip + 'en')

    df['base_cand1_ADJ.'] = base_adj
    df['base_cand2_ADJ.'] = base_adj_nouml
    df['base_cand3_VVINF'] = base_v_en

    
    return df

def get_el_bases(df_in):
    """
    Backforms the contents of column 'lemma' to the possible stems and saves candidates as new columns.
    (Corresponds to DErivBase dVN05, dNN58)

    Arg:
        df_in: dataframe containing derivation tokens in column 'lemma'
    Returns:
        df with new columns containing potential bases
    """
    df = df_in.copy()
    
    base_n = []
    base_n_nouml = []
    base_v_en = []
    base_v_en_nouml = []
    base_v_n = []
    base_v_n_nouml = []
    
    for lemma in df['lemma']:
        base_strip = lemma[:-2].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')
        base_n.append(base_strip.capitalize())
        base_n_nouml.append(base_strip.replace('ü', 'u').replace('ä', 'a').replace('ö', 'o').capitalize() if any([uml in base_strip for uml in ['ä', 'ö', 'ü']]) else '')
        base_v_en.append(base_strip + 'en')
        base_v_en_nouml.append(base_strip.replace('ü', 'u').replace('ä', 'a').replace('ö', 'o')+'en' if any([uml in base_strip for uml in ['ä', 'ö', 'ü']]) else '')
        base_v_n.append(lemma.lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö') + 'n')
        base_v_n_nouml.append(lemma.lower().replace('ü', 'u').replace('ä', 'a').replace('ö', 'o').replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')+'n' if any([uml in base_strip for uml in ['ä', 'ö', 'ü']]) else '')

    df['base_cand1_NN'] = base_n
    df['base_cand2_NN'] = base_n_nouml
    df['base_cand3_VVINF'] = base_v_en
    df['base_cand4_VVINF'] = base_v_en_nouml
    df['base_cand5_VVINF'] = base_v_n
    df['base_cand6_VVINF'] = base_v_n_nouml
    
    return df

def get_er_bases(df_in):
    """
    Backforms the contents of column 'lemma' to the possible stems and saves candidates as new columns.
    (Corresponds to DErivBase dNN05, dVN03, dNN28, dNN11)

    Arg:
        df_in: dataframe containing derivation tokens in column 'lemma'
    Returns:
        df with new columns containing potential bases
    """
    df = df_in.copy()
    
    base_v_en = []
    base_v_en_nouml = []
    base_v_rn = []
    base_v_ln = []
    base_v_ln_nouml = []
    
    base_n_asis = []
    base_n_asis_nouml = []
    base_n_noler = []
    base_n_noler_nouml = []
    base_n_ner_en_nouml = []  # Gärtner-Garten
    base_n_rer_er = [] # Abenteurer-Abenteuer
    base_n_nor_nouml = []  # Schüler-Schule
    base_n_noaner = []  # Lutheraner-Luther
    base_n_aner_a = []  # Amerikaner-Amerika
    base_n_aner_en = []  # Sizilianer-Sizilien
    
    for lemma in df['lemma']:
    
        base_strip = lemma[:-2].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')
        base_nouml = base_strip.replace('ü', 'u').replace('ä', 'a').replace('ö', 'o')
        
        if len(base_strip) <= 1:
            base_v_en.append('')
            base_v_en_nouml.append('')
            base_v_rn.append('')
            base_v_ln.append('')
            base_v_ln_nouml.append('')
            base_n_noler.append('')
            base_n_noler_nouml.append('')
            base_n_asis.append('')
            base_n_asis_nouml.append('')
            base_n_ner_en_nouml.append('')
            base_n_rer_er.append('')
            base_n_nor_nouml.append('')
            base_n_noaner.append('')
            base_n_aner_a.append('')
            base_n_aner_en.append('')
        else:
            base_v_en.append(base_strip + 'en')
            base_v_en_nouml.append(base_nouml + 'en' if any([uml in base_strip for uml in ['ä', 'ö', 'ü']]) else '')
            base_v_rn.append(base_strip + 'n' if base_strip[-1] == 'r' else '')
            base_v_ln.append(base_strip[:-1] + 'eln' if base_strip[-1] == 'l' else '')
            base_v_ln_nouml.append(base_nouml[:-1]+'eln' if (base_strip[-1] == 'l' and any([uml in base_strip for uml in ['ä', 'ö', 'ü']])) else '')
            base_n_noler.append(base_strip[:-1].capitalize() if base_strip[-1] == 'l' else '')
            base_n_noler_nouml.append(base_nouml[:-1].capitalize() if (base_strip[-1] == 'l' and any([uml in base_strip for uml in ['ä', 'ö', 'ü']])) else '')
            base_n_asis.append(base_strip.capitalize())
            base_n_asis_nouml.append(base_nouml.capitalize() if any([uml in base_strip for uml in ['ä', 'ö', 'ü']]) else '')
            base_n_ner_en_nouml.append(base_nouml[:-1].capitalize() + 'en' if base_nouml[-1] == 'n' else '')
            base_n_rer_er.append(base_strip[:-1].capitalize() + 'er' if base_strip[-1] == 'r' else '')
            base_n_nor_nouml.append(base_nouml.capitalize() + 'e')
            base_n_noaner.append(base_strip[:-2].capitalize() if base_strip[-2:] == 'an' else '')
            base_n_aner_a.append(base_strip[:-2].capitalize() + 'a' if base_strip[-2:] == 'an' else '')
            base_n_aner_en.append(base_strip[:-2].capitalize() + 'en' if base_strip[-2:] == 'an' else '')

    df['base_cand1_VVINF'] = base_v_en
    df['base_cand2_VVINF'] = base_v_en_nouml
    df['base_cand3_VVINF'] = base_v_rn
    df['base_cand4_VVINF'] = base_v_ln
    df['base_cand5_VVINF'] = base_v_ln_nouml
    
    df['base_cand6_NN'] = base_n_noler
    df['base_cand7_NN'] = base_n_noler_nouml
    df['base_cand8_NN'] = base_n_asis
    df['base_cand9_NN'] = base_n_asis_nouml
    df['base_cand10_NN'] = base_n_ner_en_nouml
    df['base_cand11_NN'] = base_n_rer_er
    df['base_cand12_NN'] = base_n_nor_nouml 
    df['base_cand13_NN'] = base_n_noaner 
    df['base_cand14_NN'] = base_n_aner_a 
    df['base_cand15_NN'] = base_n_aner_en 
    return df

def get_heit_bases(df_in):
    """
    Backforms the contents of column 'lemma' (derivations on -heit) to the possible stems and saves candidates as new columns.
    (Corresponds to DErivBase dNN59, dAN02, dAN04, dAN03.)

    Arg:
        df_in: dataframe containing heit/keit tokens in column 'lemma', output of prep_sfx_df()
    Returns:
        df with two new columns: one each for base candidate 1 and candidate 2
    """

    df = df_in.copy()

    no_heit = []
    base_n_e = []

    for lemma in df['lemma']:

        base = lemma[:-4].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')
        no_heit.append(base)

        # If the suffix was -heit, there might be a subtracted -e from the base (as in Weisheit, Trägheit, Feigheit)
        e = base + 'e' if lemma[-4] == 'h' else ''
        base_n_e.append(e)

    no_ig = [base[:-2] if base[-2:] == 'ig' else '' for base in no_heit]
    end = [base+'d' if base[-2:] == 'en' else '' for base in no_heit]

    df['base_cand1_ADJ.'] = no_heit
    df['base_cand2_ADJ.'] = no_ig
    df['base_cand3_ADJ.'] = end
    df['base_cand4_ADJ.'] = base_n_e

    return df.copy()

def get_ie_bases(df_in):
    """
    Backforms the contents of column 'lemma' to the possible stems and saves candidates as new columns.
    (Corresponds to DErivBase dAN16)

    Arg:
        df_in: dataframe containing derivation tokens in column 'lemma'
    Returns:
        df with new columns containing potential bases
    """
    df = df_in.copy()
    
    base_adj_asis = []
    base_adj_isch = []
    
    for lemma in df['lemma']:
        base_strip = lemma[:-2].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')
        base_adj_asis.append(base_strip)
        base_adj_isch.append(base_strip + 'isch')

    df['base_cand1_ADJ.'] = base_adj_asis
    df['base_cand2_ADJ.'] = base_adj_isch
    
    return df

def get_ik_bases(df_in):
    """
    Backforms the contents of column 'lemma' to the possible stems and saves candidates as new columns.
    (Corresponds to DErivBase dAN12)

    Arg:
        df_in: dataframe containing derivation tokens in column 'lemma'
    Returns:
        df with new columns containing potential bases
    """
    df = df_in.copy()
    df['base_cand1_ADJ.'] = [lemma[:-2].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö') + 'isch' for lemma in df['lemma']]
    return df

def get_iker_bases(df_in):
    """
    Backforms the contents of column 'lemma' to the possible stems and saves candidates as new columns.
    (Corresponds to DErivBase dAN13)

    Arg:
        df_in: dataframe containing derivation tokens in column 'lemma'
    Returns:
        df with new columns containing potential bases
    """
    df = df_in.copy()
    df['base_cand1_ADJ.'] = [lemma[:-4].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö') + 'isch' for lemma in df['lemma']]
    return df

def get_ikum_bases(df_in):
    """
    Backforms the contents of column 'lemma' to the possible stems and saves candidates as new columns.
    (Corresponds to DErivBase dAN14)

    Arg:
        df_in: dataframe containing derivation tokens in column 'lemma'
    Returns:
        df with new columns containing potential bases
    """
    df = df_in.copy()
    df['base_cand1_ADJ.'] = [lemma[:-4].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö') + 'isch' for lemma in df['lemma']]
    return df

def get_ismus_bases(df_in):
    """
    Backforms the contents of column 'lemma' to the possible stems and saves candidates as new columns.
    (Corresponds to DErivBase dAN09)

    Arg:
        df_in: dataframe containing derivation tokens in column 'lemma'
    Returns:
        df with new columns containing potential bases
    """
    df = df_in.copy()
    
    base_adj_asis = []
    base_adj_ell = []
    base_adj_ar = []
    base_adj_isch = []
    base_adj_zisch = []
    
    base_n_asis = []
    base_n_ik = []
    
    for lemma in df['lemma']:
        base_strip = lemma[:-5].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')

        base_adj_asis.append(base_strip)
        base_adj_ell.append(base_strip[:-2] + 'ell' if base_strip[-2:] == 'al' else '')
        base_adj_ar.append(base_strip[:-2] + 'är' if base_strip[-2:] == 'ar' else '')
        base_adj_isch.append(base_strip + 'isch')
        base_adj_zisch.append(base_strip[:-2] + 'isch' if base_strip[-2:] == 'iz' else '')
        
        base_n_asis.append(base_strip.capitalize())
        base_n_ik.append(base_strip.capitalize() + 'ik')
        
    df['base_cand1_ADJ.'] = base_adj_asis
    df['base_cand2_ADJ.'] = base_adj_ell
    df['base_cand3_ADJ.'] = base_adj_ar
    df['base_cand4_ADJ.'] = base_adj_isch
    df['base_cand5_ADJ.'] = base_adj_zisch
    
    df['base_cand6_NN'] = base_n_asis
    df['base_cand6_NN'] = base_n_ik
    
    return df

def get_ist_bases(df_in):
    """
    Backforms the contents of column 'lemma' to the possible stems and saves candidates as new columns.
    (Corresponds to DErivBase dAN15, dVN31, dNN08)

    Arg:
        df_in: dataframe containing derivation tokens in column 'lemma'
    Returns:
        df with new columns containing potential bases
    """
    df = df_in.copy()
    
    base_adj_asis = []
    base_adj_isch = []
    base_adj_add_isch = []
    base_v_ieren = []
    base_n_asis = []
    base_n_misc = []
    
    lastvowel = {'Poliz':'Polizei', 'Pian':'Piano', 'Gitarr':'Gitarre', 'Anarch':'Anarchie', 'Propagand':'Propaganda'}
    
    for lemma in df['lemma']:
        base = lemma.lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')
        base_strip = lemma[:-3].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')

        base_adj_asis.append(base_strip)
        base_adj_isch.append(base_strip + 'isch')
        base_adj_add_isch.append(base + 'isch')
        base_v_ieren.append(base_strip + 'ieren')
        
        n = base_strip.capitalize()
        base_n_asis.append(n)
        
        base_n_misc.append( lastvowel[n] if n in lastvowel.keys() else '' )
            
        
    df['base_cand1_ADJ.'] = base_adj_asis
    df['base_cand2_ADJ.'] = base_adj_isch
    df['base_cand3_ADJ.'] = base_adj_add_isch
    df['base_cand4_VVINF'] = base_v_ieren
    df['base_cand5_NN'] = base_n_misc
    df['base_cand6_NN'] = base_n_asis
    
    return df
    
def get_itaet_bases(df_in):
    """
    Backforms the contents of column 'lemma' to the possible stems and saves candidates as new columns.
    (Corresponds to DErivBase dAN16)

    Arg:
        df_in: dataframe containing derivation tokens in column 'lemma'
    Returns:
        df with new columns containing potential bases
    """
    df = df_in.copy()
    
    base_adj_asis = []
    base_adj_ell = []
    base_adj_zisch = []
    base_adj_bel = []
    
    for lemma in df['lemma']:
      # Even though -ität only has four characters, I guess ä is represented in Py2 with two chars?
      # So we have to remove the last five in order to get the correct base.
        base_strip = lemma[:-5].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')

        base_adj_asis.append(base_strip)
        base_adj_ell.append(base_strip[:-2] + 'ell' if base_strip[-2:] == 'al' else '')
        base_adj_zisch.append(base_strip[:-2] + 'isch' if base_strip[-2:] == 'iz' else '')
        base_adj_bel.append(base_strip[:-3] + 'bel' if base_strip[-3:] == 'bil' else '')

    df['base_cand1_ADJ.'] = base_adj_asis
    df['base_cand2_ADJ.'] = base_adj_ell
    df['base_cand3_ADJ.'] = base_adj_bel
    df['base_cand4_ADJ.'] = base_adj_zisch
    
    return df

def get_ition_bases(df_in):
    """
    Backforms the contents of column 'lemma' to the possible stems and saves candidates as new columns.
    (Corresponds to DErivBase dVN34)

    Arg:
        df_in: dataframe containing derivation tokens in column 'lemma'
    Returns:
        df with new columns containing potential bases
    """
    df = df_in.copy()
    
    base_v_ieren = []
    base_v_nieren = []
    
    for lemma in df['lemma']:
        base_strip = lemma[:-5].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')

        base_ieren = base_strip + 'ieren'
        base_nieren = base_strip[:-1] + 'nieren' if base_strip[-1] == 's' else ''
        
        base_v_ieren.append(base_ieren)
        base_v_nieren.append(base_nieren)

    df['base_cand1_VVINF'] = base_v_ieren
    df['base_cand2_VVINF'] = base_v_nieren
    
    return df

def get_ium_bases(df_in):
    """
    Backforms the contents of column 'lemma'.

    Arg:
        df_in: dataframe containing derivations in column 'lemma'
    Returns:
        df with new columns containing potential bases
    """
    
    df = df_in.copy()
    df['base_cand1_VVINF'] = [lemma[:-3].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')+'ieren' for lemma in df['lemma']]
    df['base_cand2_NN'] = [lemma[:-3] for lemma in df['lemma']]
    return df

def get_ling_bases(df_in):
    """
    Backforms the contents of column 'lemma' to the possible stems and saves candidates as new columns.
    (Corresponds to DErivBase dAN08, extended)

    Arg:
        df_in: dataframe containing derivation tokens in column 'lemma'
    Returns:
        df with new columns containing potential bases
    """
    df = df_in.copy()
    
    base_adj_asis = []
    base_adj_nouml = []
    base_adj_e = []
    base_v_en = []
    base_v_en_nouml = []
    base_v_n = []
    
    for lemma in df['lemma']:
        base_strip = lemma[:-4].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')

        adj_asis = base_strip
        adj_nouml = base_strip.replace('ü', 'u').replace('ä', 'a').replace('ö', 'o') if any([uml in base_strip for uml in ['ä', 'ö', 'ü']]) else ''
        adj_e = base_strip + 'e' if base_strip in {'feig', 'träg', 'weis'} else ''
        
        v_en = adj_asis + 'en'
        v_en_nouml = adj_nouml + 'en' if adj_nouml != '' else ''
        v_n = adj_asis + 'n' if base_strip[-1] in ['r', 'l'] else ''
        
        base_adj_asis.append(adj_asis)
        base_adj_nouml.append(adj_nouml)
        base_adj_e.append(adj_e)
        base_v_en.append(v_en)
        base_v_n.append(v_n)
        base_v_en_nouml.append(v_en_nouml)

    df['base_cand1_ADJ.'] = base_adj_asis
    df['base_cand2_ADJ.'] = base_adj_nouml
    df['base_cand3_ADJ.'] = base_adj_e
    df['base_cand4_VVINF'] = base_v_en
    df['base_cand5_VVINF'] = base_v_n
    df['base_cand6_VVINF'] = base_v_en_nouml
    
    return df

def get_nis_bases(df_in):
    """
    Backforms the contents of column 'lemma' (derivations on -nis) to the possible stems and saves candidates as new columns.
    The base might be adjectives (geheim/Geheimnis, finster/Finsternis, gefangen/Gefängnis) or verbs (ärgern/Ärgernis, erleben/Erlebnis).
    (Corresponds to DErivBase dAN05, dAN07, dVN06)

    Arg:
        df_in: dataframe containing -nis tokens in column 'lemma'
    Returns:
        df with five new columns containing potential bases
    """

    df = df_in.copy()

    # Drop all lines with words that we know we don't want.
    donotwant = {'Tennis', 'Anis', 'Penis', 'Dennis'}
    df = df[~df['lemma'].isin(donotwant)]
    
    base_adj_simp = []
    base_adj_ppart = []
    base_v_n = []
    base_v_en = []
    base_v_nen = []

    # Irregular stems and where blanket removal of ä is invalid
    irreg = {'kennt':'kenn', 'gedacht':'gedenk', 'arger':'ärger', 'bedrang':'bedräng', 'stand':'steh', 'verhang':'verhäng', 'versaum':'versäum'}

    for lemma in df['lemma']:

        # Strip off -nis
        base = lemma[:-3].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')

        # Get rid of N1 in hyphenated compounds (this is now done in the orig query script, can remove for later samples)
        if len(base.split('-')) > 1:
            base = base.split('-')[1]

        # De-umlautify ä
        base = base.replace('ä', 'a')

        # Check if the base is in the dict of irregular stems; if it is, replace it with the value string.
        for key, val in irreg.items():
            if key in base:
                base = base.replace(key, val)

        # If the base is a simplex adjective, it should now be formed.
        base_adj_simp.append(base)

        # If it's a past participle like gefangen, add -en if it starts with 'ge'
        ppart = base+'en' if base[:2] == 'ge' else ''
        base_adj_ppart.append(ppart)

        # If it's a verb, the infinitive could be formed with -en (the simplest case), or also -n or -nen.
        n = base+'n' if base[-2:] == 'er' else ''
        nen = base+'nen' if base[-1] == 'g' or base[-2:] == 'ch' else ''
        base_v_n.append(n)
        base_v_nen.append(nen)
        base_v_en.append(base+'en')

    # Add columns to df.
    df['base_cand1_ADJ.'] = base_adj_simp
    df['base_cand2_ADJ.'] = base_adj_ppart
    df['base_cand3_VVINF'] = base_v_n
    df['base_cand4_VVINF'] = base_v_en
    df['base_cand5_VVINF'] = base_v_nen

    return df

def get_schaft_bases(df_in):
    """
    Backforms the contents of column 'lemma' (derivations on -schaft) to the possible stems and saves candidates as new columns.
    The bases can be verbs  (pflegen/Pflegschaft, wandern/Wanderschaft), nouns (Anwalt/Anwaltschaft, Friede/Friedschaft), or adjs (schwanger/Schwangerschaft).
    (Corresponds to DErivBase dNN04, dAN06, dVN11)

    Arg:
        df_in: dataframe containing -schaft tokens in column 'lemma'
    Returns:
        df with seven new columns containing potential bases
    """

    df = df_in.copy()

    base_v_en = []
    base_v_ern = []
    base_n_e = []
    base_n_simp = []
    base_n_no_n = []
    base_n_no_en = []
    base_adj = []

    # Irregular/umlauty stems
    irreg = {'Brüder':'Bruder'}

    for lemma in df['lemma']:

        base = lemma[:-6]

        # Check if the base is in the dict of irregular stems; if it is, replace it with the value string.
        for key, val in irreg.items():
            if key in base:
                base = base.replace(key, val)

        # For the nouns, just strip the suffix off and add -e
        base_n_simp.append(base)
        base_n_e.append(base + 'e')
        no_n = base[:-1] if base[-2:] == 'en' else ''
        no_en = base[:-2] if base[-2:] == 'en' else ''
        base_n_no_n.append(no_n)
        base_n_no_en.append(no_en)

        # For verbs and adjectives, lowercase.
        base = base.lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')
        base_adj.append(base)

        # For verbs, add 'en', and if base ends in 'r', then 'n'
        base_v_en.append(base + 'en')
        try:
            rn = base + 'n' if base[-1] == 'r' else ''
        except: # if index error because base is too short
            rn = ''
        base_v_ern.append(rn)

    df['base_cand1_NN'] = base_n_simp
    df['base_cand2_NN'] = base_n_e
    df['base_cand3_ADJ.'] = base_adj
    df['base_cand4_VVINF'] = base_v_en
    df['base_cand5_VVINF'] = base_v_ern
    df['base_cand6_NN'] = base_n_no_n
    df['base_cand7_NN'] = base_n_no_en

    return df.copy()

def get_pfx_bases(df_in, pfx_string):
    """
    Backforms the contents of column 'lemma' (derivations with a prefix) to the possible stems and saves candidates as new columns.
    The bases are all nouns. 
    NB: This function is NOT for Ge-; it has its own.
    
    Args:
        df_in: dataframe containing -schaft tokens in column 'lemma'
        pfx_string: string representing the current prefix w hyphen, e.g., 'Um-'
    Returns:
        df with seven new columns containing potential bases
    """

    df = df_in.copy()
    pfx_len = len(pfx_string)-1
    df['base_cand1_NN'] = [lemma[pfx_len:].capitalize() for lemma in df['lemma']]
    return df

def get_ung_bases(df_in):
    """
    Backforms the contents of column 'lemma' to the possible stems and saves candidates as new columns.
    (Corresponds to DErivBase dVN07)

    Arg:
        df_in: dataframe containing derivation tokens in column 'lemma'
    Returns:
        df with new columns containing potential bases
    """
    df = df_in.copy()
    
    base_v_en = []
    base_v_rn = []
    base_v_ln = []

    for lemma in df['lemma']:
        base = lemma[:-3].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')
        base_v_en.append(base + 'en')
        base_v_rn.append(base + 'n' if base[-1] == 'r' else '')
        base_v_ln.append(base[:-1] + 'eln' if base[-1] == 'l' else '')
    
    df['base_cand1_VVINF'] = base_v_en
    df['base_cand2_VVINF'] = base_v_rn
    df['base_cand3_VVINF'] = base_v_ln
    
    return df
    
def get_bases_w_ieren(df_in, sfx_str):
    """
    Backforms the contents of column 'lemma' to verbs that end in -ieren and saves candidates as new columns.
    Used for DErivBase dVN21 (-and), dVN23 (-ant), dVN25 (-age), dVN26 (-anz), dVN29 (-ent), dVN30 (-enz), dVN28 (-ement), 
    dVN27 (-atur), dVN33 (-itur), dVN37 (-iment), dVN38 (-eur), dVN39 (-ation)

    Arg:
        df_in: dataframe containing derivations in column 'lemma'
        sfx_str: string containing the suffix, e.g., '-and'
    Returns:
        df with new columns containing potential bases
    """
    
    df = df_in.copy()
    len_sfx = len(sfx_str)-1
    df['base_cand1_VVINF'] = [lemma[:-len_sfx].lower().replace('Ü', 'ü').replace('Ä', 'ä').replace('Ö', 'ö')+'ieren' for lemma in df['lemma']]
    
    return df

def get_unique_base_cands(df):
    """
    Identifies the unique bases from a df containing potential base candidates. Returns them, along with their POS, in a new df.

    Arg:
        df: pandas df, output of get_X_bases()
    Returns:
        Pandas df with three columns: lemma, unique_candidates, and pos.
    """

    # Convert the dataframe into a format that will make deduplication easier while maintaining lemma-base mapping.
    cands = pd.melt(df.drop(['morph', 'doc.id', 'doc.url', 's.idx', 'word'], axis=1),
                    id_vars = 'lemma',
                    var_name = 'base_type',
                    value_name = 'unique_candidates')

    cands['pos'] = cands['base_type'].str.split('_').str[-1]
    cands.drop(['base_type'], axis=1, inplace=True)

    # Drop duplicates in 'cand' and 'pos', get rid of rows containing only '' as candidate, and return.
    cands = cands.drop_duplicates(subset=['unique_candidates', 'pos']).reset_index(drop=True)
    cands = cands[cands['unique_candidates'] != '']

    return cands

def add_cql(df):
    """
    Creates CQL queries for use in SeaCOW for each base candidate and its respective POS

    Arg:
        df: pandas df, output of get_unique_base_cands()
    Returns:
        Pandas df
    """

    # Define a function that takes a pd.Series (a row of df) and formats the two data points as CQL query string.
    def cql_row(row):
        return '[lemma="{0}" & tag="{1}"] within <s/>'.format(row['unique_candidates'], row['pos'])

    # Apply that function to df and return.
    df['cql'] = df[['unique_candidates', 'pos']].apply(cql_row, axis=1)
    return df

def get_cql_from_bases(bases_df):
    """
    Uses other methods in this class to identify unique candidate bases and format them as CQL queries, ready for SeaCOW.

    Arg:
        bases_df: dataframe containing original data along with the generated bases for the given data and sfx
    Returns:
        A list of strings representing the CQL queries
    """
    unique_cands_df = get_unique_base_cands(bases_df)
    cql_list = add_cql(unique_cands_df)
    return cql_list

def get_bases(df_raw, sfx):
    """
    Frontend function that calls the get_X_bases() method of whatever suffix is passed in and returns
    whatever was returned.

    Arg:
        df_raw: dataframe containing derivations in column 'lemma', read in from file.
        sfx: string representing the suffix (must be in predetermined list of suffixes)
    Returns:
        pandas df containing original data along with the generated bases for the given data and sfx.
    """
    IEREN_SFXS = {'-age', '-and', '-ant', '-anz', '-ation', '-atur', '-ement', '-end', '-ent', '-enz', '-eur', '-iment', '-iteur', '-itur'}
    INDIV_SFXS = {'-ament', '-ateur', '-ator', '-e', '-el', '-er', '-heit', '-ie', '-ik', '-iker', '-ikum', '-ismus', '-ist', '-itaet', '-ition', '-ium', '-ling', '-nis', '-schaft', '-ung'}
    
    # Check that input suffix is valid.
    if sfx not in set.union(IEREN_SFXS, INDIV_SFXS):
        raise ValueError('Invalid suffix, please enter one of the following forms: %s' % str(set.union(IEREN_SFXS, INDIV_SFXS)))
        
    # Clean up the data using prep_sfx_df().
    df = prep_sfx_df(df_raw)
    	
    if sfx in IEREN_SFXS:
        bases = get_bases_w_ieren(df, sfx)
    elif sfx in {'-ament', '-ateur'}:
        bases = get_ament_ateur_bases(df)
    elif sfx == '-ator':
        bases = get_ator_bases(df)
    elif sfx == '-e':
        bases = get_e_bases(df)
    elif sfx == '-el':
        bases = get_el_bases(df)
    elif sfx == '-er':
        bases = get_er_bases(df)
    elif sfx == '-heit':
        bases = get_heit_bases(df)
    elif sfx == '-ie':
        bases = get_ie_bases(df)
    elif sfx == '-ik':
        bases = get_ik_bases(df)
    elif sfx == '-iker':
        bases = get_iker_bases(df)
    elif sfx == '-ikum':
        bases = get_ikum_bases(df)
    elif sfx == '-ismus':
        bases = get_ismus_bases(df)
    elif sfx == '-ist':
        bases = get_ist_bases(df)
    elif sfx == '-itaet':
        bases = get_itaet_bases(df)
    elif sfx == '-ition':
        bases = get_ition_bases(df)
    elif sfx == '-ium':
        bases = get_ium_bases(df)
    elif sfx == '-ling':
        bases = get_ling_bases(df)
    elif sfx == '-nis':
        bases = get_nis_bases(df)
    elif sfx == '-schaft':
        bases = get_schaft_bases(df)
    elif sfx == '-ung':
        bases = get_ung_bases(df)

    return bases

