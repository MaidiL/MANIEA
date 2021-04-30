import sys
import time
from numpy import *
import feather
import csv
import networkx as nx
import matplotlib.pyplot as plt
import random

fo = feather.read_dataframe('edd_singh.otu_table.clean.feather')
title = fo.columns

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

def dataDevisionOtuData(Path1,Path2):
    print("Data being divided into case and control...")
    fm = feather.read_dataframe(Path1)
    fo = feather.read_dataframe(Path2)
    healthIndex = []
    eddIndex = []
    healthMetadata = fm
    eddMetadata = fm
    healthOtudata = fo
    eddOtudata = fo
    count = 0
    state = fm.DiseaseState
    for item in state:
        if item == 'H':
            healthIndex.append(count)
        elif item == 'EDD':
            eddIndex.append(count)
        count += 1

    print("Dividing datasets...\nPlease wait...\n")

    # Data to be discarded by each data group
    healthDropIndex = []
    eddDropIndex = []

    for index in healthIndex:
        eddDropIndex.append(index)

    for index in eddIndex:
        healthDropIndex.append(index)

    eddDropIndex.sort()
    healthDropIndex.sort()

    healthDropLen = len(healthDropIndex)
    eddDropLen = len(eddDropIndex)

    # Traverse from back to front
    for i in range(healthDropLen -1, -1, -1):
        dropIndex = healthDropIndex[i]
        healthOtudata = healthOtudata.drop(index = [dropIndex])
        healthMetadata = healthMetadata.drop(index = [dropIndex])

    for i in range(eddDropLen -1, -1, -1):
        dropIndex = eddDropIndex[i]
        eddOtudata = eddOtudata.drop(index = [dropIndex])
        eddMetadata = eddMetadata.drop(index = [dropIndex])

    return healthMetadata,eddMetadata,healthOtudata,eddOtudata

# Preprocess the data
def loadDataSet(data):
    print("Start loading dataset...\n")
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
                if item[column] != 0:
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

dedge = []
uedge = []

with open("Singh-Fpgrowth-UndirectedEdges-0.65-0.75.csv",'r',newline='') as f:
    f_csv = csv.reader(f)
    for item in f_csv:
        uedge.append(item)

with open("Singh-Fpgrowth-DirectedEdges-0.65-0.75.csv",'r',newline='') as f:
    f_csv = csv.reader(f)
    for item in f_csv:
        dedge.append(item)

# Undirected graph
U = nx.Graph()
# Directed graph
D = nx.DiGraph()

U.add_edges_from(uedge)
D.add_edges_from(dedge)

# Composite Chart
G = nx.compose(U, D)
position = nx.kamada_kawai_layout(G)
print('position', position)

# set the size of the nodes
de = G.degree()
de2 = []
for item in de:
    de2.append(item[1] * 50)

print("G.nodes", G.nodes())
print("G.degree", G.degree())
print("G.edges", len(G.edges), G.edges())

family = {}
for node in G.nodes():
    familyNm = findFamily(node,title)
    familyNm = familyNm.lstrip('f__')
    if familyNm in family.keys():
        family[familyNm].append(node)
    else:
        family[familyNm]=[node]
print('family',family)
print('family.keys',family.keys())

def randomcolor():
    colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
    color = ""
    for i in range(6):
        color += colorArr[random.randint(0,14)]
    return "#"+color

random.seed(0)

colorDict = {}
for item in family.keys():
    colorDict[item] = randomcolor()

print('colorDict',colorDict)

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
nx.draw(G, with_labels=True, node_size=de2, node_color=colorLs1, alpha=1, linewidths=3,
        font_color='black', edge_color='grey', width=1, style='dashdot',
        pos=position, font_size=10)
nx.draw(D, with_labels=False, node_size=de2, node_color=colorLs2, alpha=1, linewidths=3,
        font_color='black', edge_color='cornflowerblue', width=1, style='solid',
        pos=position)
plt.savefig('Singh-Fpgrowth-NetworkModel-DrawNetworks.png',dpi=500)
plt.show()

print("0. Overall network indexes")
print("Number of nodes", len(nx.degree(G)))
print("Number of Undirected edges", U.size())
print("Number of positive directed edges", D.size())

l = G.size()

print("1. Density")
print(l * 2 / (len(nx.degree(G)) * (len(nx.degree(G)) - 1)))

print("2. Average connectivity")
degreeSum = 0
for item in nx.degree(G):
    degreeSum += item[1]
aveConnectivity = degreeSum / len(nx.degree(G))
print(aveConnectivity)

print("3. Average length of the shortest paths")
print(nx.average_shortest_path_length(G))

print("4. Average clustering coefficient")
cluserSum = 0
for item in nx.clustering(G).values():
    cluserSum += item
aveClustering = cluserSum / len(nx.degree(G))
print(aveClustering)

print("5. Transitivity")
print(nx.transitivity(G))

print("Individual network indexes")
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

healthMetadata, eddMetadata, healthOtudata, eddOtudata = \
    dataDevisionOtuData('edd_singh.metadata.clean.feather', 'edd_singh.otu_table.clean.feather')
dataset1, title1 = loadDataSet(healthOtudata)

abundanceLs = []

nodeLs=[]
for item in G.nodes():
    for i in range(len(title1)):
        if item in title1[i]:
            nodeLs.append(i)
print('nodeLs',nodeLs)

for node in nodeLs:
    nodeNm = 0
    for item in dataset1:
        if node in item:
            nodeNm += 1
    abundanceLs.append(nodeNm)

print("abundanceLs",abundanceLs)
print("Average occurrence frequency",mean(abundanceLs))

# Save the network indexes for individual nodes
with open('Singh-Fpgrowth-IndexforIndividualNodes.csv', 'w', newline='') as f:
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
for i in range(10):
    rankDict[ls[i][0]] += 10-i

ls = sorted(nx.betweenness_centrality(G).items(), key=itemgetter(1), reverse=True)
for i in range(10):
    rankDict[ls[i][0]] += 10-i

ls = sorted(nx.closeness_centrality(G).items(), key=itemgetter(1), reverse=True)
for i in range(10):
    rankDict[ls[i][0]] += 10-i

ls = sorted(nx.eigenvector_centrality(G).items(), key=itemgetter(1), reverse=True)
for i in range(10):
    rankDict[ls[i][0]] += 10-i

ls = sorted(nx.clustering(G).items(), key=itemgetter(1), reverse=True)
for i in range(10):
    rankDict[ls[i][0]] += 10-i

print('rankDict',sorted(rankDict.items(), key=itemgetter(1), reverse=True))