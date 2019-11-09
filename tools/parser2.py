import mygene
import pickle
import os
import gc 

biogrid_save = "biogrid.pkl"

if(os.path.exists(biogrid_save)):
    interact = pickle.load(open(biogrid_save,'rb'))
else:
    interact = {}
    with open("BIOGRID-ORGANISM-Homo_sapiens-3.5.177.tab2.txt",'r') as f:
        next(f)
        it = 0
        for line in f:
            cols = line.strip().split("\t")
            
            officialA, officialB = ( cols[7], cols[8] )
            synonymA, synonymB = ( cols[9], cols[10] )
            
            it+=1
            print(it)

            if officialA not in interact:
                interact[officialA] = [officialB]
                for synB in synonymB.split("|"):
                    interact[officialA].append(synB)
            else:
                interact[officialA].append(officialB)
                for synB in synonymB.split("|"):
                    interact[officialA].append(synB)

            for synA in synonymA.split("|"):
                if(synA not in interact):
                    interact[synA] = [officialB]
                    for synB in synonymB.split("|"):
                        interact[synA].append(synB)
                else:
                    interact[synA].append(officialB)
                    for synB in synonymB.split("|"):
                        interact[officialA].append(synB)

    pickle.dump( interact , open(biogrid_save, "wb" ) , pickle.HIGHEST_PROTOCOL)