import mygene
import pickle
import os
import gc


mg = mygene.MyGeneInfo()


save1 = "save1.pkl"
save2 = "save2.pkl"
save3 = "save3.pkl"


saveQueryA = "BiogridQueryA.pkl"
saveQueryB = "BiogridQueryB.pkl"


#create with HI-union
if(os.path.exists(save1)):
    ensg_gene, interaction = pickle.load(open(save1,'rb'))
else:
    ensg_gene = {}
    interaction = {}

    with open("HI-union.tsv",'r') as f:
        for line in f:
            cols = line.strip().split("\t")
            p1 = cols[0]
            p2 = cols[1]
            
            if(p1 not in ensg_gene):
                infos = mg.getgenes(p1, fields="uniprot")
                try:
                    genes = []
                    infos = infos[0]["uniprot"]
                    for k,v in infos.items():
                        if(not isinstance(v,list)):
                            genes.append(v)
                        else:
                            genes.extend(v)
                    ensg_gene[p1] = genes

                except KeyError:
                    print(infos)
                    pass
            
            # orienté p1 -> p2 , pas dans les 2 sens 
            if(p1 not in interaction):
                interaction[p1] = [p2]
            else:
                interaction[p1].append(p2)
        
    pickle.dump( (ensg_gene,interaction) , open(save1, "wb" ) )



#enhance with HuRI dataset
if(os.path.exists(save2)):
    ensg_gene, interaction = pickle.load(open(save2,'rb'))
else:
    with open("HuRI.tsv",'r') as f:
        for line in f:
            cols = line.strip().split("\t")
            p1 = cols[0]
            p2 = cols[1]
            
            if(p1 not in ensg_gene):
                infos = mg.getgenes(p1, fields="uniprot")
                try:
                    genes = []
                    infos = infos[0]["uniprot"]
                    for k,v in infos.items():
                        if(not isinstance(v,list)):
                            genes.append(v)
                        else:
                            genes.extend(v)
                    ensg_gene[p1] = genes

                except KeyError:
                    print(infos)
                    pass
            
            if(p1 not in interaction):
                interaction[p1] = [p2]
            else:
                interaction[p1].append(p2)
        
    pickle.dump( (ensg_gene,interaction) , open(save2, "wb" ) )


#enhance with BIOGRID dataset
if(os.path.exists(saveQueryA) and os.path.exists(saveQueryB)):
    p1_infos = pickle.load(open(saveQueryA,'rb'))
    p2_infos = pickle.load(open(saveQueryB,'rb'))
else:
    PA = []
    PB = []
    with open("BIOGRID-ORGANISM-Homo_sapiens-3.5.177.tab2.txt",'r') as f:
        next(f)
        for line in f:
            cols = line.strip().split("\t")
            p1 = cols[1]
            p2 = cols[2]

            PA.append(p1)
            PB.append(p2)

    p1_infos = mg.getgenes(PA, fields="entrezgene, uniprot, ensembl.gene")
    PA = None
    gc.collect()
    pickle.dump( p1_infos , open(saveQueryA, "wb" ) , pickle.HIGHEST_PROTOCOL)

    p2_infos = mg.getgenes(PB, fields="entrezgene, uniprot, ensembl.gene")
    PB = None
    gc.collect()
    pickle.dump( p2_infos , open(saveQueryB, "wb" ) , pickle.HIGHEST_PROTOCOL)


if(os.path.exists(save3)):
    ensg_gene, interaction = pickle.load(open(save3,'rb'))
else:

    for infos in zip(p1_infos, p2_infos):
        try:
            if(isinstance(infos[0]["ensembl"],list)):
                pA = [ g["gene"] for g in infos[0]["ensembl"] ]
            else:
                pA = [ infos[0]["ensembl"]["gene"] ]

            if(isinstance(infos[1]["ensembl"],list)):
                pB = [ g["gene"] for g in infos[1]["ensembl"] ]
            else:
                pB = [ infos[1]["ensembl"]["gene"] ]
        except KeyError: 
            print(infos)
            continue
            
        for p1 in pA:
            if(p1 not in ensg_gene):
                try:
                    genes = []
                    infos = infos[0]["uniprot"]
                    for k,v in infos.items():
                        if(not isinstance(v,list)):
                            genes.append(v)
                        else:
                            genes.extend(v)

                    if( p1 not in ensg_gene ):
                        ensg_gene[p1] = genes
                    else:
                        #extend without duplicates
                        mylist = ensg_gene[p1].copy().extend(genes)
                        mylist = list(dict.fromkeys(mylist))
                        ensg_gene[p1] = mylist

                except KeyError:
                    print(infos)
                    pass

            for p2 in pB:
                if(p1 not in interaction):
                    interaction[p1] = [p2]
                else:
                    interaction[p1].append(p2)

        pickle.dump( (ensg_gene,interaction) , open(save3, "wb" ) )
