import os
import pickle
import pandas as pd
import numpy as np

TFs = [ "H2AZ", "H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me2", "H3K4me3","H3K79me2","H3K9ac","H3K9me1",
        "H3K9me3","H4K20me1","CREBBP","CBX2","CBX8","CHD1","CHD7","SETDB1","EZH2","HDAC1","HDAC2","HDAC6","CBX3",
        "KDM5C","KDM1A","WHSC1","PHF8","KDM5B","POL2b","RBBP5","RNAPIIS5P","SUZ12"]

save = "biogrid.pkl"

if(os.path.exists(save)):
    interaction = pickle.load(open(save,'rb'))
else:
    exit(1)



mat = np.zeros( (len(TFs), len(TFs)) )
notfound = []
found = 0 
for i in range(len(TFs)):
    pA = TFs[i]
    ni = 0
    if pA not in interaction:
        notfound.append(pA)
        continue
    for pB in interaction[pA]:
        if pB in TFs:
            found += 1 
            j = TFs.index(pB)
            mat[i,j] += 1
            ni += 1

    if(ni != 0):
        mat[i]




print(mat)
print(notfound)
print(found)

with open("experimental_graph.txt",'w') as f:
    autoregulator = 0 
    all_inter = 0
    for i in range(len(TFs)):
        for j in range(len(TFs)):
            c = mat[i,j]
            if(c != 0):
                #f.write(TFs[i] + '  ' + TFs[j] + ' ' +str(1) +'\n')
                f.write(str(i) + '  ' + str(j) +'\n')

                all_inter += 1
                if(i==j):
                    autoregulator += 1
    
    print(autoregulator, all_inter)


    average_inter = all_inter / (len(TFs)-len(notfound))
    estimate_edges = average_inter * len(TFs)
    
    print(estimate_edges)




