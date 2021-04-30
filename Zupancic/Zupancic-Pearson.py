import sys
import time
from numpy import *
import feather
import csv
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import random
import numpy as np

type = sys.getfilesystemencoding()

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

    print("healthIndex:\n", healthIndex)
    print("obIndex:\n", obIndex)
    print("nonIndex:\n", nonIndex)

    healthLen = len(healthIndex)
    obLen = len(obIndex)
    nonLen = len(nonIndex)

    print("healthLen:\n", healthLen)
    print("obLen:\n", obLen)
    print("nonLen:\n", nonLen)
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

    print('nonDropIndex',nonDropIndex)
    print('obDropIndex',obDropIndex)
    print('healthDropIndex',healthDropIndex)

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
def loadDataset3(data):
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
            valueLs = []
            for column in genusRecord[genus]:
                valueLs.append(item[column])
                newValue[countrow][countcolumn] = np.max(valueLs)
            countrow += 1

    newTitle.insert(0, 'ID')
    value = newValue
    title = newTitle
    print(len(value[0]), len(title))
    print("Dataset has been combined into Genus level...\n")

    for i in range(len(value)):
        for j in range(len(value[0])):
            if value[i][j] == 0:
                value[i][j] = 0.1

    return value, title

def getFamily(dataset, title):
    allTitle = dataset.columns
    allTitleList = []
    familyRecord = {}
    for item in allTitle[1:]:
        ls = item.split(";")
        allTitleList.append(ls)
    for genus in title:
        for item in allTitleList:
            for name in item:
                if genus == name[3:]:
                    index = item.index(name)
                    family = item[index - 1]
                    family = family.lstrip('f__')
                    if family not in familyRecord.keys():
                        familyRecord[family] = [genus]
                    else:
                        if genus not in familyRecord[family]:
                            familyRecord[family].append(genus)
    return familyRecord

if __name__=='__main__':
    # Dividing dataset
    healthMetadata, obMetadata, nonMetadata, healthOtudata, obOtudata, nonOtudata = \
        dataDevisionOtuData('ob_zupancic.metadata.clean.feather', 'ob_zupancic.otu_table.clean.feather')

    # Set minimum Pearson threshold
    minCor = 0.37

    # Preprocessing dataset
    dataset, title = loadDataset3(healthOtudata)
    fo = feather.read_dataframe('ob_zupancic.otu_table.clean.feather')
    title1 = fo.columns

    for item in dataset:
        del item[0]
    del title[0]
    print("dataset",len(dataset),dataset)
    print("title",len(title),title)
    abundanceVector = []
    genusNum = len(title)
    for i in range(genusNum):
        abundanceVector.append([])
        for item in dataset:
            abundanceVector[-1].append(item[i])
    print('abundanceVector',len(abundanceVector),abundanceVector)

    # calculating correlation coefficient value
    corRecorder = []
    time_start = time.time()
    for number1 in range(len(title)):
        print("progress",number1)
        for number2 in range(number1+1,len(title)):
            x = pd.Series(abundanceVector[number1])
            y = pd.Series(abundanceVector[number2])
            corValue = x.corr(y, method='pearson')
            corRecorder.append([title[number1],title[number2],corValue])
    time_end = time.time()
    print("corRecorder",len(corRecorder),corRecorder)

    # recording correlations satisfying the threshold
    posCorRecorder = []
    negCorRecorder = []
    for item in corRecorder:
        if item[2] > minCor:
            posCorRecorder.append(item)
        if item[2] < -1 * minCor:
            negCorRecorder.append(item)

    # constructing network model
    posuedge = []
    neguedge = []
    time_start2 = time.time()
    for item in posCorRecorder:
        posuedge.append([item[0].lstrip('g__'), item[1].lstrip('g__')])
    for item in negCorRecorder:
        neguedge.append([item[0].lstrip('g__'), item[1].lstrip('g__')])

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

    # setting layout of nodes
    position = nx.spring_layout(U, k=0.3)

    def findFamily(gName, bioTitle):
        for bio in bioTitle:
            if 'g__' + gName + ';' in bio:
                for charNm in range(len(bio) - 3):
                    if bio[charNm:charNm + 3] == 'f__':
                        charNm1 = charNm
                    if bio[charNm:charNm + 3] == 'g__':
                        charNm2 = charNm
                fName = bio[charNm1:charNm2 - 1]
                break
        return fName

    family = {}
    for node in U.nodes():
        familyNm = findFamily(node, title1)
        familyNm = familyNm.lstrip('f__')
        if familyNm in family.keys():
            family[familyNm].append(node)
        else:
            family[familyNm] = [node]

    # setting color of nodes according to the genus the microorganisms belong to
    def randomcolor():
        colorArr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F']
        color = ""
        for i in range(6):
            color += colorArr[random.randint(0, 14)]
        return "#" + color

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

    nx.draw(U1, with_labels=False, node_size=de2, node_color=colorLs1, alpha=1, linewidths=3,
            font_color='black', edge_color='cornflowerblue', width=1, style='solid',
            pos=position, font_size=11)

    nx.draw(U2, with_labels=False, node_size=de2, node_color=colorLs2, alpha=1, linewidths=3,
            font_color='black', edge_color='red', width=1, style='solid',
            pos=position)

    plt.savefig('Zupancic-Pearson-NetworkModel-' + str(minCor) + '.png', dpi=500)
    plt.show()
    time_end2 = time.time()

    with open('Zupancic-Pearson-AssociationRules-'+str(minCor)+'.csv', 'w', newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerow(["Correlation calculation time", time_end - time_start])
        f_csv.writerow(["Network construction time", time_end2 - time_start2])
        f_csv.writerow(["A", "B", "Cor"])
        f_csv.writerow(["posCorRecorder"])
        f_csv.writerows(posCorRecorder)
        f_csv.writerow(["negCorRecorder"])
        f_csv.writerows(negCorRecorder)

    with open('Zupancic-Pearson-PositiveEdges-' + str(minCor) + '.csv', 'w', newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerows(posuedge)
    with open('Zupancic-Pearson-NegativeEdges-' + str(minCor) + '.csv', 'w', newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerows(neguedge)

    print("Correlation calcualtion time", time_end - time_start)
    print("Network construction time", time_end2 - time_start2)