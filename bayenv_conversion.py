import pandas as pd
import os

files = ['TER_angsd_polysites.mafs', 'LOM_angsd_polysites.mafs']
data = []

for csv in files:
    print(csv)
    frame = pd.read_csv(csv, sep="\t")
    frame['filename'] = csv.split('_')[0]
    data.append(frame)
    
def init_proc_file(data):
    #data = data.loc[data['nInd'] == 20] #only keep sites that has data from every individual
    #print(data.shape[0])
    data = data.rename(columns={"chromo": "nuc", "position": "major", "major": "minor", "minor": "alg", "knownEM": "al_freq", "pK-EM": "pval"})
    data.index.name = 'chr'
    data.reset_index(inplace=True)
    data['nuc'] = data['nuc'].apply(str)
    data['chr'] = data['chr'].apply(str)
    print(type(data['chr'][1]))
    data['pos'] = data[['chr', 'nuc']].apply(lambda x: ''.join(x), axis=1) #create position so that we can merge dfs later
    df = data[["pos","al_freq", "filename"]]
    print(df.shape[1])
    df["al_count"] = df["al_freq"] * 20 #Note that allele counts are not integers. This is because of the uncertainty retained by ANGSD. But here, we have to get rid of that. So, let's round to closest integer.
    df["al_count"] = df.al_count.round().astype(int)
    df["al_count2"] = 20 - df["al_count"] 
    df = df[["pos","filename", "al_count", "al_count2"]]
    pop = df["filename"][2]
    print(pop)
    df = df.rename(columns={"al_count": "al_count"+pop, "al_count2": "al_count2"+pop})
    df = df.drop(['filename'], axis=1)
    return df
    
ter_df = init_proc_file(data[0])
lom_df = init_proc_file(data[1])

pd.merge(ter_df, lom_df, on="pos")
