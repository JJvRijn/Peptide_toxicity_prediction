"""
This is a boilerplate pipeline 'feature_calculation'
generated using Kedro 0.18.7
"""
import pandas as pd
from itertools import product
#grabs all the unique sequences so that features are calculated just once
def unique_peptides (complete_data) -> list:
    list_of_peptides = complete_data.seq
    return list(set(list_of_peptides))

#split AA in Peptide
def split_peptide(peptide):
    return [AA for AA in peptide]

def AAC_ss(seq):
    NAA = ["G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P",
           "V", "I", "C", "Y", "H", "R", "N", "D", "T"]
    AAC = []
    split_pep = split_peptide(seq)
    for x in range(len(NAA)):
        AACount = split_pep.count(NAA[x])
        AANumber = len(seq)
        AAC.append(AACount/AANumber)
    return AAC

#Amino Acid Composition Dataframe
def AAC_df(seqs):
    colnames = ["%G", "%A", "%L", "%M", "%F", "%W", "%K", "%Q", "%E", "%S", "%P",
                "%V", "%I", "%C", "%Y", "%H", "%R", "%N", "%D", "%T"]
    AAC_df = pd.DataFrame(columns = colnames)
    for y in range(len(seqs)):
        AAC_pep = AAC_ss(seqs[y])
        AAC_df.loc[y] = AAC_pep
    AAC_df["seq"] = seqs
    return AAC_df

#di Amino Acid Composition single sequence
def DAAC_ss(seq, i=1):
    test2 = split_peptide(seq)
    NAA = ["G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P",
           "V", "I", "C", "Y", "H", "R", "N", "D", "T"]
    z = i + 1      #2 means next neighbor, 3 skips 1, 4 skips 2 etc
    all_comb_2 = list(product(NAA, NAA))
    count = [0]*len(all_comb_2)
    comp_DAAC = [0]*len(all_comb_2)
    for j in range(len(all_comb_2)):
        for v in range(len(test2) - 1):
            if list(all_comb_2[j]) == test2[v:v+z]:
                count[j] = count[j] + 1
    for r in range(len(all_comb_2)):
        comp_DAAC[r] = count[r]/sum(count)
    return comp_DAAC

#di Amino Acid Composition for dataframe
def DAAC_df(seqs, i=1):
    NAA = ["G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P",
           "V", "I", "C", "Y", "H", "R", "N", "D", "T"]
    all_comb_2 = list(map(str, list(product(NAA, NAA))))
    colnames = colnames = ["%" + s[2] + s[7] for s in all_comb_2]
    DAAC_df = pd.DataFrame(columns = colnames)
    for y in range(len(seqs)):
        if len(seqs[y]) > 2:
            comp_DAAC = DAAC_ss(seqs[y], i)
            DAAC_df.loc[y] = comp_DAAC
        else:
            DAAC_df.loc[y] = [0] * len(colnames)
    DAAC_df["seq"] = seqs
    return DAAC_df

#multi Amino Acid Composition single sequence
def MAAC_ss(seq, n=3, i=1):
    test2 = split_peptide(seq)
    NAA = ["G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P",
           "V", "I", "C", "Y", "H", "R", "N", "D", "T"]
    all_comb_n = list(product(NAA, repeat=n))
    count = [0]*len(all_comb_n)
    comp_MAAC = [0]*len(all_comb_n)
    #checks if Dipeptide, cuz than seperation works
    if i != 1:
        if n == 2:
            i = i
        else:
            i = 1
            print("Higher order only works for dipeptides. \nComputed with i =1.")
    z = i + n- 1 #2 means next neighbor, 3 skips 1, 4 skips 2 etc
    for j in range(len(all_comb_n)):
        for v in range(len(test2) - 1):
            if list(all_comb_n[j]) == test2[v:v+z]:
                count[j] = count[j] + 1
    for r in range(len(all_comb_n)):
        if sum(count) != 0:
            comp_MAAC[r] = count[r]/sum(count)
    return comp_MAAC

def MAAC_df(seqs, n=3, i=1):
    NAA = ["G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P",
           "V", "I", "C", "Y", "H", "R", "N", "D", "T"]
    all_comb_n = list(map(str, list(product(NAA, repeat = n))))
    colnames = ["%" + s[2] for s in all_comb_n]
    for o in range(n-1):
        holder = [s[2+(5*(o+1))] for s in all_comb_n]
        for a in range(len(holder)):
            colnames[a] = colnames[a] + holder[a]
    MAAC_df = pd.DataFrame(columns = colnames)
    for y in range(len(seqs)):
        if len(seqs[i]) > 3:
            comp_MAAC = MAAC_ss(seqs[y], n=n,  i=i)
            MAAC_df.loc[y] = comp_MAAC
        else:
            MAAC_df.loc[y] = [0] * len(colnames)
    MAAC_df["seq"] = seqs
    return MAAC_df

def compact_fp_ss(seq):
    phi = ["V", "I", "L", "F", "W", "Y", "M", "A"] # = Hydrophobic
    omega = ["F", "W", "Y", "H"] # = Aromatic
    psi = ["V", "I", "L", "M"] # = Aliphatic (non-aromatic)
    pi = ["P", "G", "A", "S"] # = small
    zeta = ["S", "T", "H", "N", "Q", "E", "D", "K", "R"] # = Hydrophilic
    s_groups = ["C", "U"]
    pos = ["K", "R", "H"] # = Positively charged
    neg = ["D", "E"]
    cat = [phi, omega, psi, pi, zeta, pos, neg, s_groups]
    comp_fp = []
    split_pep = split_peptide(seq)
    for x in range(len(split_pep)):
        for y in range(len(cat)):
            comp_fp.append(float(split_pep[x] in cat[y]))
    return comp_fp

#compact fingerprint for dataframe
def compact_fp_df(seqs):
    cat_names = ["phi", "omega", "psi", "pi", "zeta", "pos", "neg", "s_group"]
    #length = seqs.str.len()
    max_length =len( max(seqs, key=len))
    colnames = []
    for a in range(max_length):
        colname = ([z + str(a) for z in cat_names])
        colnames = colnames + colname
    df_comp = pd.DataFrame(columns = colnames)
    for b in range(len(seqs)):
        fp_of_pep = compact_fp_ss(seqs[b])
        fp_of_pep = fp_of_pep + [0]*((max_length*8) - len(fp_of_pep))
        df_comp.loc[b] = fp_of_pep
    df_comp["seq"] = seqs
    return df_comp


def AA_matrix(seq):
    import numpy as np
    NAA = ["G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P",
           "V", "I", "C", "Y", "H", "R", "N", "D", "T"]
    split_pep = split_peptide(seq)
    AA_matrix = np.zeros((50, 20))
    for j in range(len(split_pep)):
        for o in range(len(NAA)):
            AA_matrix[j][o] = int((NAA[o] == split_pep[j]))
    return AA_matrix

def AA_matrix_df(seqs):
    import pandas as pd
    df_comp = pd.DataFrame(columns=["AA_matrix"])
    for n in range(len(seqs)):
        df_comp.loc[n] = [AA_matrix(seqs[n])]
    df_comp["seq"] = seqs
    return df_comp

def BLOSUM62(seq):
    import numpy as np
    split_pep = split_peptide(seq)
    BLOSUM62_matrix = np.zeros((len(seq), 20))
    blosum_values = {
        "A" : [4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0],
        "R" : [-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3],
        "N" : [-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3],
        "D" : [-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3],
        "C" : [0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1 ],
        "Q" : [-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2],
        "E" : [-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2],
        "G" : [0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3],
        "H" : [-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3],
        "I" : [-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3],
        "L" : [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1],
        "K" : [-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2],
        "M" : [-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1],
        "F" : [-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1],
        "P" : [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2],
        "S" : [1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2],
        "T" : [0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0],
        "W" : [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3],
        "Y" : [-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1],
        "V" : [0, -3, -3, -3, -1 ,-2 ,-2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4 ],
    }
    for i in range(len(split_pep)):
        BLOSUM62_matrix[i] = blosum_values.get(split_pep[i])
    return BLOSUM62_matrix

def BLOSUM62_df(seqs):
    import pandas as pd
    df_comp = pd.DataFrame(columns = ["Blosum62"])
    for v in range(len(seqs)):
        df_comp.loc[v] = [BLOSUM62(seqs[v])]
    df_comp["seq"] = seqs
    return df_comp


def chunk(in_string,num_chunks):
    chunk_size = len(in_string)//num_chunks
    if len(in_string) % num_chunks: chunk_size += 1
    iterator = iter(in_string)
    for _ in range(num_chunks):
        accumulator = list()
        for _ in range(chunk_size):
            try: accumulator.append(next(iterator))
            except StopIteration: break
        yield ''.join(accumulator)
        
def distrubution(seq):
    NAA = ["G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P",
           "V", "I", "C", "Y", "H", "R", "N", "D", "T"]
    dist = []
    chunks = list(chunk(seq, 4))
    for a in chunks:
        split_pep = split_peptide(a)
        for n in range(len(NAA)):
            AACount = split_pep.count(NAA[n])
            AANumber = len(seq)
            dist.append(AACount/AANumber)
    return dist

def distrubution_df(seqs):
    colnames = ["25%_G%", "25%_A%", "25%_L%", "25%_M%", "25%_F%", "25%_W%", "25%_K%", "25%_Q%",
                "25%_E%", "25%_S%", "25%_P%", "25%_V%", "25%_I%", "25%_C%", "25%_Y%", "25%_H%",
                "25%_R%", "25%_N%", "25%_D%", "25%_T%", "50%_G%", "50%_A%", "50%_L%", "50%_M%",
                "50%_F%", "50%_W%", "50%_K%", "50%_Q%", "50%_E%", "50%_S%", "50%_P%", "50%_V%",
                "50%_I%", "50%_C%", "50%_Y%", "50%_H%", "50%_R%", "50%_N%", "50%_D%", "50%_T%",
                "75%_G%", "75%_A%", "75%_L%", "75%_M%", "75%_F%", "75%_W%", "75%_K%", "75%_Q%",
                "75%_E%", "75%_S%", "75%_P%", "75%_V%", "75%_I%", "75%_C%", "75%_Y%", "75%_H%",
                "75%_R%", "75%_N%", "75%_D%", "75%_T%", "100%_G%", "100%_A%", "100%_L%",
                "100%_M%", "100%_F%", "100%_W%", "100%_K%", "100%_Q%", "100%_E%", "100%_S%",
                "100%_P%", "100%_V%", "100%_I%", "100%_C%", "100%_Y%", "100%_H%", "100%_R%",
                "100%_N%", "100%_D%", "100%_T%"]
    df_comp = pd.DataFrame(columns = colnames)
    for a in range(len(seqs)):
        df_comp.loc[a] = distrubution(seqs[a])
    df_comp["seq"] = seqs
    return df_comp