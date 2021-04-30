from numpy import *
import feather
import csv
import networkx as nx
import matplotlib.pyplot as plt

fo = feather.read_dataframe('crc_baxter.otu_table.clean.feather')
title1 = fo.columns

def findFamily(gName, bioTitle):
    for bio in bioTitle:
        if 'g__'+gName+';' in bio:
            for charNm in range(len(bio)-3):
                if bio[charNm:charNm+3] == 'f__':
                    charNm1 = charNm
                if bio[charNm:charNm+3] == 'g__':
                    charNm2 = charNm
            fName = bio[charNm1:charNm2-1]
            break
    return fName

# Divide the data by health status of samples
def dataDevisionOtuData(Path1,Path2):
    print("Data being divided into case and control...")
    fm = feather.read_dataframe(Path1)
    fo = feather.read_dataframe(Path2)
    healthIndex = []
    crcIndex = []
    noncrcIndex = []
    healthMetadata = fm
    crcMetadata = fm
    noncrcMetadata = fm
    healthOtudata = fo
    crcOtudata = fo
    noncrcOtudata = fo
    count = 0
    state = fm.DiseaseState
    for item in state:
        if item == 'H':
            healthIndex.append(count)
        elif item == 'CRC':
            crcIndex.append(count)
        elif item == 'nonCRC':
            noncrcIndex.append(count)
        count += 1

    print("Dividing datasets...\nPlease wait...\n")

    # Data to be discarded by each data group
    healthDropIndex = []
    noncrcDropIndex = []
    crcDropIndex = []

    for index in healthIndex:
        noncrcDropIndex.append(index)
        crcDropIndex.append(index)

    for index in crcIndex:
        healthDropIndex.append(index)
        noncrcDropIndex.append(index)

    for index in noncrcIndex:
        healthDropIndex.append(index)
        crcDropIndex.append(index)

    noncrcDropIndex.sort()
    crcDropIndex.sort()
    healthDropIndex.sort()

    healthDropLen = len(healthDropIndex)
    crcDropLen = len(crcDropIndex)
    noncrcDropLen = len(noncrcDropIndex)

    # Traverse from back to front
    for i in range(healthDropLen -1, -1, -1):
        dropIndex = healthDropIndex[i]
        healthOtudata = healthOtudata.drop(index = [dropIndex])
        healthMetadata = healthMetadata.drop(index = [dropIndex])

    for i in range(crcDropLen -1, -1, -1):
        dropIndex = crcDropIndex[i]
        crcOtudata = crcOtudata.drop(index = [dropIndex])
        crcMetadata = crcMetadata.drop(index = [dropIndex])

    for i in range(noncrcDropLen -1, -1, -1):
        dropIndex = noncrcDropIndex[i]
        noncrcOtudata = noncrcOtudata.drop(index = [dropIndex])
        noncrcMetadata = noncrcMetadata.drop(index=[dropIndex])

    return healthMetadata,crcMetadata,noncrcMetadata,healthOtudata,crcOtudata,noncrcOtudata

# Preprocess the data
def loadDataSet(data):
    print("Start preprocessing dataset...\n")
    f = data
    # Title is the header
    title = f.columns
    # Value is a sparse matrix
    value = f.values
    recordDic = {}
    # Record which columns of data each genus contains
    genusRecord = {}
    count = 0
    for item in title:
        ls = item.split(";")
        for word in ls:
            if "g__" in word:
                if word in genusRecord.keys():
                    genusRecord[word].append(count)
                else:
                    genusRecord[word] = []
                    genusRecord[word].append(count)
        count += 1
    del genusRecord['g__']
    # Title list after merging into genus
    newTitle = []
    # Value matrix after merging into genus
    newValue = []

    for key in genusRecord.keys():
        newTitle.append(key)

    for i in range(len(value)):
        newValue.append([value[i][0]])
        for j in range(len(newTitle)):
            newValue[i].append(0)

    # Merge to genus level
    countcolumn = 0
    for genus in newTitle:
        countrow = 0
        countcolumn += 1
        for item in value:
            for column in genusRecord[genus]:
                if item[column] >= 5:
                    newValue[countrow][countcolumn] = 1
                    break
            countrow += 1

    newTitle.insert(0,'index')
    value = newValue
    title = newTitle

    # Eliminated data with all 1 or all 0, because the existence of these microorganisms is little affected by other microorganisms
    uselessColumns = []
    usefulColumns = []
    maximum = len(value[0])
    for j in range(1, maximum):
        if value[0][j] == 0:
            for i in range(len(value)):
                if value[i][j] != 0:
                    usefulColumns.append(j)
                    break
        if value[0][j] != 0:
            for i in range(len(value)):
                if value[i][j] == 0:
                    usefulColumns.append(j)
                    break
    for j in range(1,maximum):
        if j not in usefulColumns:
            uselessColumns.append(j)
    uselessColumns.sort()
    uselessLen = len(uselessColumns)

    for i in range(uselessLen-1, -1, -1):
        del title[uselessColumns[i]]
        for item in value:
            del item[uselessColumns[i]]

    for item in value:
        recordDic[item[0]] = []
        i = 0
        maximum = len(item)
        for i in range(maximum):
            if item[i] == 1:
                recordDic[item[0]].append(i)
                continue
    # Create data list
    data = []
    for keys in recordDic:
        data.append(recordDic[keys])
    return data,title

title = ['index', 'g__Clostridium_XI', 'g__Faecalibacterium', 'g__Prevotella', 'g__Ruminococcus', 'g__Roseburia', 'g__Akkermansia', 'g__Lachnoanaerobaculum', 'g__Alistipes', 'g__Lachnospiracea_incertae_sedis', 'g__Phascolarctobacterium', 'g__Acidaminococcus', 'g__Haemophilus', 'g__Dialister', 'g__Escherichia/Shigella', 'g__Ruminococcus2', 'g__Coprococcus', 'g__Anaerostipes', 'g__Eggerthella', 'g__Clostridium_XlVa', 'g__Parabacteroides', 'g__Oscillibacter', 'g__Clostridium_IV', 'g__Allisonella', 'g__Bifidobacterium', 'g__Clostridium_sensu_stricto', 'g__Odoribacter', 'g__Gemmiger', 'g__Streptococcus', 'g__Flavonifractor', 'g__Clostridium_XlVb', 'g__Barnesiella', 'g__Collinsella', 'g__Anaerofilum', 'g__Anaerofustis', 'g__Megasphaera', 'g__Methanobrevibacter', 'g__Lactobacillus', 'g__Clostridium_XVIII', 'g__Dorea', 'g__Paraprevotella', 'g__Veillonella', 'g__Victivallis', 'g__Hespellia', 'g__Peptoniphilus', 'g__Sutterella', 'g__Eubacterium', 'g__Butyricicoccus', 'g__Pseudoflavonifractor', 'g__Fusobacterium', 'g__Subdoligranulum', 'g__Erysipelotrichaceae_incertae_sedis', 'g__Peptostreptococcus', 'g__Parasutterella', 'g__Rothia', 'g__Methanosphaera', 'g__Desulfovibrio', 'g__Methylobacterium', 'g__Sporobacter', 'g__Hydrogenoanaerobacterium', 'g__Enterococcus', 'g__Anaerosporobacter', 'g__Gracilibacter', 'g__Acetanaerobacterium', 'g__Lachnobacterium', 'g__Oxalobacter', 'g__Clostridium_III', 'g__Paenibacillus', 'g__Parvimonas', 'g__Gemella', 'g__Acetivibrio', 'g__Leuconostoc', 'g__Cerasicoccus', 'g__Turicibacter', 'g__Actinomyces', 'g__Neisseria', 'g__Anaerococcus', 'g__Pseudomonas', 'g__Enterobacter', 'g__Defluviitalea', 'g__Anaerotruncus', 'g__Olsenella', 'g__Oribacterium', 'g__Weissella', 'g__Methanomassiliicoccus', 'g__Porphyromonas', 'g__Bilophila', 'g__Butyricimonas', 'g__Pyramidobacter', 'g__Mogibacterium', 'g__Listeria', 'g__Coprobacillus', 'g__Cloacibacillus', 'g__Slackia', 'g__Catenibacterium', 'g__Ralstonia', 'g__Megamonas', 'g__Hallella', 'g__Succinivibrio', 'g__Corynebacterium', 'g__Alloscardovia', 'g__Granulicatella', 'g__Vampirovibrio', 'g__Alloprevotella', 'g__Robinsoniella', 'g__Enterorhabdus', 'g__Lactonifactor', 'g__Gardnerella', 'g__Finegoldia', 'g__Marvinbryantia', 'g__Peptococcus', 'g__Acetitomaculum', 'g__Holdemania', 'g__Acetobacterium', 'g__Howardella', 'g__Anaeroplasma', 'g__Citrobacter', 'g__Campylobacter', 'g__Actinobacillus', 'g__Lactococcus', 'g__Anaerophaga', 'g__Solobacterium', 'g__Deinococcus', 'g__Succiniclasticum', 'g__Streptophyta', 'g__Gordonibacter', 'g__Acinetobacter', 'g__Anaerovorax', 'g__Eggerthia', 'g__Moryella', 'g__Ethanoligenens', 'g__Oxobacter', 'g__Anaerorhabdus', 'g__Selenomonas', 'g__Papillibacter', 'g__Clostridium_XII', 'g__Dietzia', 'g__Atopobium', 'g__Bacillus', 'g__Saccharibacteria_genera_incertae_sedis', 'g__Sarcina', 'g__Salmonella', 'g__Flavobacterium', 'g__Abiotrophia', 'g__Allobaculum', 'g__Leptotrichia', 'g__Pantoea', 'g__Staphylococcus', 'g__Pediococcus', 'g__Varibaculum', 'g__Vagococcus', 'g__Pasteurella']

uedge = []
dedge = []

with open("Baxter-Apriori-UndirectedEdges-0.8-0.9.csv",'r',newline='') as f:
    f_csv = csv.reader(f)
    for item in f_csv:
        uedge.append(item)

with open("Baxter-Apriori-DirectedEdges-0.8-0.9.csv",'r',newline='') as f:
    f_csv = csv.reader(f)
    for item in f_csv:
        dedge.append(item)

# Undirected graph
U = nx.Graph()
# Directed graph
D = nx.DiGraph()

U.add_edges_from(uedge)
D.add_edges_from(dedge)

print("Number of undirected edges", len(U.edges), U.edges())
print("Number of directed edges", len(dedge))

# Composite Chart
G = nx.compose(U, D)
position = nx.kamada_kawai_layout(G)
print('position', position)

nodeLs = []
for node in G.nodes():
    nodeLs.append(node)
for node in D.nodes():
    if node in nodeLs:
        continue
    else:
        nodeLs.append(node)

de = G.degree()
de2 = []
for item in de:
    de2.append(item[1] * 25)
de3 = []
for item in de:
    de3.append(item[1] * 40)

print("G.nodes", G.nodes())
print("G.degree", G.degree())
print("G.edges", len(G.edges), G.edges())

family = {}
for node in G.nodes():
    familyNm = findFamily(node,title1)
    familyNm = familyNm.lstrip('f__')
    if familyNm in family.keys():
        family[familyNm].append(node)
    else:
        family[familyNm]=[node]
print('family',family)
print('family.keys',family.keys())

colorDict = {'Ruminococcaceae': '#E7DF71', 'Rikenellaceae': '#5987FD', 'Lachnospiraceae': '#E586AF', 'Bacteroidaceae': '#F49353', 'Erysipelotrichaceae': '#D2AD5F'}
colorLs1 = []
for node in G.nodes:
    for item in family.keys():
        if node in family[item]:
            colorLs1.append(colorDict[item])
            break
colorLs2 = []
for node in D.nodes:
    for item in family.keys():
        if node in family[item]:
            colorLs2.append(colorDict[item])
            break

position = nx.kamada_kawai_layout(G)
print("position", position)
plt.figure(figsize=(7.5, 6))
plt.tight_layout()
nx.draw(G, with_labels=True, node_size=de3, node_color=colorLs1, alpha=1, linewidths=3,
        font_color='black', edge_color='grey', width=1, style='dashdot',
        pos=position, font_size=11)
nx.draw(D, with_labels=False, node_size=de3, node_color=colorLs2, alpha=1, linewidths=3,
        font_color='black', edge_color='cornflowerblue', width=1, style='solid',
        pos=position)
plt.savefig('Baxter-Apriori-NetworkModel-DrawNetworks.png',dpi=500)
plt.show()

print("0. Overall network indexes")
print("Number of nodes", len(nx.degree(G)))
print("Number of Undirected edges", U.size())
print("Number of positive directed edges", D.size())

l = G.size()

print("1. Density")
print(l*2 / (len(nx.degree(G))*(len(nx.degree(G))-1)))

print("2. Average connectivity")
degreeSum = 0
for item in nx.degree(G):
    degreeSum += item[1]
aveConnectivity = degreeSum/len(nx.degree(G))
print(aveConnectivity)

print("3. Average length of the shortest paths")
print(nx.average_shortest_path_length(G))

print("4. Average clustering coefficient")
cluserSum = 0
for item in nx.clustering(G).values():
    cluserSum += item
aveClustering = cluserSum/len(nx.degree(G))
print(aveClustering)

print("5. Transitivity")
print(nx.transitivity(G))

print("Network indexes for individual nodes")
from operator import itemgetter
print("1. Degree")
degreeDict = {}
for item in G.degree():
    degreeDict[item[0]] = 0
for item in U.degree():
    degreeDict[item[0]] += item [1]
for item in D.degree():
    degreeDict[item[0]] += item [1]
print(sorted(degreeDict.items(),key=itemgetter(1),reverse=True))

print("2. Betweenness")
print(sorted(nx.betweenness_centrality(G).items(),key=itemgetter(1),reverse=True))

print("3. Closeness centrality")
print(sorted(nx.closeness_centrality(G).items(),key=itemgetter(1),reverse=True))

print("4. Eigenvector centrality")
print(sorted(nx.eigenvector_centrality(G).items(),key=itemgetter(1),reverse=True))

print("5. Clustering coefficient")
print(sorted(nx.clustering(G).items(),key=itemgetter(1),reverse=True))

print("6. Average occurrence frequency")

healthMetadata, crcMetadata, noncrcMetadata, healthOtudata, crcOtudata, noncrcOtudata = \
    dataDevisionOtuData('crc_baxter.metadata.clean.feather', 'crc_baxter.otu_table.clean.feather')
dataset1, title1 = loadDataSet(healthOtudata)

abundanceLs = []

nodeLs=[]
for item in G.nodes():
    for i in range(len(title1)):
        if item in title1[i]:
            nodeLs.append(i)

for node in nodeLs:
    nodeNm = 0
    for item in dataset1:
        if node in item:
            nodeNm += 1
    abundanceLs.append(nodeNm)

print("Average occurrence frequency",mean(abundanceLs))

# Save the network indexes for individual nodes
with open('Baxter-Apriori-IndexforIndividualNodes.csv', 'w', newline='') as f:
    f_csv = csv.writer(f)
    f_csv.writerow(['degree'])
    for item in sorted(degreeDict.items(), key=itemgetter(1), reverse=True):
        f_csv.writerow([item[0],item[1]])
    f_csv.writerow(['betweenness'])
    for item in sorted(nx.betweenness_centrality(G).items(), key=itemgetter(1), reverse=True):
        f_csv.writerow([item[0],item[1]])
    f_csv.writerow(['Closeness centrality'])
    for item in sorted(nx.closeness_centrality(G).items(),key=itemgetter(1),reverse=True):
        f_csv.writerow([item[0], item[1]])
    f_csv.writerow(['Eigenvector centrality'])
    for item in sorted(nx.eigenvector_centrality(G).items(),key=itemgetter(1),reverse=True):
        f_csv.writerow([item[0], item[1]])
    f_csv.writerow(['Clustering coefficient'])
    for item in sorted(nx.clustering(G).items(),key=itemgetter(1),reverse=True):
        f_csv.writerow([item[0], item[1]])

rankDict = {}
for node in G.nodes:
    rankDict[node] = 0

ls = sorted(degreeDict.items(),key=itemgetter(1),reverse=True)
for i in range(5):
    rankDict[ls[i][0]] += 5-i

ls = sorted(nx.betweenness_centrality(G).items(), key=itemgetter(1), reverse=True)
for i in range(5):
    rankDict[ls[i][0]] += 5-i

ls = sorted(nx.closeness_centrality(G).items(), key=itemgetter(1), reverse=True)
for i in range(5):
    rankDict[ls[i][0]] += 5-i

ls = sorted(nx.eigenvector_centrality(G).items(), key=itemgetter(1), reverse=True)
for i in range(5):
    rankDict[ls[i][0]] += 5-i

ls = sorted(nx.clustering(G).items(), key=itemgetter(1), reverse=True)
for i in range(5):
    rankDict[ls[i][0]] += 5-i

print('rankDict',sorted(rankDict.items(), key=itemgetter(1), reverse=True))