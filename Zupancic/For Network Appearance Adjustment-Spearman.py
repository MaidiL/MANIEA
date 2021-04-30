import sys
import time
from numpy import *
import feather
import csv
import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations
import pandas as pd
import random

random.seed(0)

fo = feather.read_dataframe('ob_zupancic.otu_table.clean.feather')
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
    obIndex = []
    nonIndex = []
    healthMetadata = fm
    obMetadata = fm
    nonMetadata = fm
    healthOtudata = fo
    obOtudata = fo
    nonOtudata = fo
    count = 0
    state = fm.DiseaseState
    for item in state:
        if item == 'H':
            healthIndex.append(count)
        elif item == 'OB':
            obIndex.append(count)
        else:
            nonIndex.append(count)
        count += 1

    print("Dividing datasets...\nPlease wait...\n")

    # Data to be discarded by each data group
    healthDropIndex = []
    nonDropIndex = []
    obDropIndex = []

    for index in healthIndex:
        nonDropIndex.append(index)
        obDropIndex.append(index)

    for index in obIndex:
        healthDropIndex.append(index)
        nonDropIndex.append(index)

    for index in nonIndex:
        healthDropIndex.append(index)
        obDropIndex.append(index)

    nonDropIndex.sort()
    obDropIndex.sort()
    healthDropIndex.sort()

    healthDropLen = len(healthDropIndex)
    obDropLen = len(obDropIndex)
    nonDropLen = len(nonDropIndex)

    # Traverse from back to front
    for i in range(healthDropLen -1, -1, -1):
        dropIndex = healthDropIndex[i]
        healthOtudata = healthOtudata.drop(index = [dropIndex])
        healthMetadata = healthMetadata.drop(index = [dropIndex])

    for i in range(obDropLen -1, -1, -1):
        dropIndex = obDropIndex[i]
        obOtudata = obOtudata.drop(index = [dropIndex])
        obMetadata = obMetadata.drop(index = [dropIndex])

    for i in range(nonDropLen -1, -1, -1):
        dropIndex = nonDropIndex[i]
        nonOtudata = nonOtudata.drop(index = [dropIndex])
        nonMetadata = nonMetadata.drop(index=[dropIndex])

    return healthMetadata,obMetadata,nonMetadata,healthOtudata,obOtudata,nonOtudata

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

posuedge = []
neguedge = []

with open("Zupancic-Spearman-PositiveEdges-0.39.csv",'r',newline='') as f:
    f_csv = csv.reader(f)
    for item in f_csv:
        posuedge.append(item)

with open("Zupancic-Spearman-NegativeEdges-0.39.csv",'r',newline='') as f:
    f_csv = csv.reader(f)
    for item in f_csv:
        neguedge.append(item)

print("posuedge",len(posuedge),posuedge)
print("neguedge",len(neguedge),neguedge)

# Undirected graph1 for positive edges
U1 = nx.Graph()
U1.add_edges_from(posuedge)

# Undirected graph2 for negative edges
U2 = nx.Graph()
U2.add_edges_from(neguedge)

# Composite Chart
U = nx.compose(U1, U2)

de = U.degree()
de2 = []
for item in de:
    de2.append(item[1] * 50)

family = {}
for node in U.nodes():
    familyNm = findFamily(node,title1)
    familyNm = familyNm.lstrip('f__')
    if familyNm in family.keys():
        family[familyNm].append(node)
    else:
        family[familyNm]=[node]

def randomcolor():
    colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
    color = ""
    for i in range(6):
        color += colorArr[random.randint(0,14)]
    return "#"+color

colorDict = {}
for item in family.keys():
    colorDict[item] = randomcolor()

colorLs1 = []
for node in U1.nodes:
    for item in family.keys():
        if node in family[item]:
            colorLs1.append(colorDict[item])
            break

colorLs2 = []
for node in U2.nodes:
    for item in family.keys():
        if node in family[item]:
            colorLs2.append(colorDict[item])
            break

position = nx.spring_layout(U, k=0.3)
print('position', position)

nx.draw(U1, with_labels=False, node_size=de2, node_color=colorLs1, alpha=1, linewidths=3,
        font_color='black', edge_color='cornflowerblue', width=1, style='solid',
        pos=position, font_size=11)

nx.draw(U2, with_labels=False, node_size=de2, node_color=colorLs2, alpha=1, linewidths=3,
        font_color='black', edge_color='orangered', width=1, style='solid',
        pos=position)

plt.savefig('Zupancic-Spearman-NetworkModel-DrawNetworks.png',dpi=500)
plt.show()

print("Overall network indexes")
print("0. Number of nodes and edges")
print("Number of nodes",len(nx.degree(U)))
print("Number of undirected positive edges",U1.size())
print("Number of undirected negative edges",U2.size())

l = U.size()

print("1. Density")
print(l*2 / (len(nx.degree(U))*(len(nx.degree(U))-1)))

print("2.Average connectivity")
degreeSum = 0
for item in nx.degree(U):
    degreeSum += item[1]
aveConnectivity = degreeSum/len(nx.degree(U))
print(aveConnectivity)

print("3. Average length of the shortest paths")
# print(nx.average_shortest_path_length(U))

print("4. Average clustering coefficient")
cluserSum = 0
for item in nx.clustering(U).values():
    cluserSum += item
aveClustering = cluserSum/len(nx.degree(U))
print(aveClustering)

print("5. Transitivity")
print(nx.transitivity(U))

print("6. Average occurrence frequency")
healthMetadata, obMetadata, nonMetadata, healthOtudata, obOtudata, nonOtudata = \
        dataDevisionOtuData('ob_zupancic.metadata.clean.feather', 'ob_zupancic.otu_table.clean.feather')
dataset1, title1 = loadDataSet(healthOtudata)
abundanceLs = []

nodeLs=[]
for item in U.nodes():
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