import os
import pickle
import pandas as pd
import numpy as np

save1 = "save1.pkl"
save2 = "save2.pkl"
save3 = "save3.pkl"

save = save3

if(os.path.exists(save)):
    ensg_gene, interaction = pickle.load(open(save,'rb'))
else:
    exit(1)



print(len(ensg_gene)  ,len(interaction))

TFs = [ "H2AZ", "H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me2", "H3K4me3","H3K79me2","H3K9ac","H3K9me1",
        "H3K9me3","H4K20me1","CREBBP","CBX2","CBX8","CHD1","CHD7","SETDB1","EZH2","HDAC1","HDAC2","HDAC6","CBX3",
        "KDM5C","KDM1A","WHSC1","PHF8","KDM5B","POL2b","RBBP5","RNAPIIS5P","SUZ12"]


TF_to_ensg = [[]] * len(TFs)

for i in range(len(TFs)):
    tf = TFs[i]
    for k,v in ensg_gene.items():
        if(tf in v):
            TF_to_ensg[i].append(k)


# aucun match HuRI & HI-union

print(TF_to_ensg)
#test si []


mat = np.zeros( (len(TFs), len(TFs)) )

for i in range(len(TFs)):
    pA = TF_to_ensg[i]
    if(len(pA) == 0):
        continue
    ni = 0
    for p1 in pA:
        for j in range(len(TFs)):
            pB = TF_to_ensg[j]
            for p2 in pB:
                if( p2 in interaction[p1] ):
                    mat[i,j] += 1
                    ni += 1

    if(ni != 0):
        mat[i] /= ni 



df = pd.read_table("../data/norm_preProcessed_32nodes.txt", sep="\t")

dico_gene = {}

for gene in df.columns:
    dico_gene[gene] = {"mean": df[gene].mean(), 'ecart-type': df[gene].var()}

# print(dico_gene)