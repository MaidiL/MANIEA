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

type = sys.getfilesystemencoding()

def eclat(prefix, items, min_support, freq_items):
    while items:
        # get itemset and corresponding transaction number from items one by one
        key, item = items.pop()
        # The number of transactions is the support of the itemset
        key_support = len(item)
        # Determine if itemsets are frequent
        if key_support >= min_support:
            # If the itemset is frequent, add the itemset and its support to the frequent itemset set freq_items
            freq_items[frozenset(sorted(prefix+[key]))] = key_support
            suffix = []
            for other_key, other_item in items:
                new_item = item & other_item
                # Determine if the candidate itemset occur frequently
                # If frequent, add it to suffix
                if len(new_item) >= min_support:
                    suffix.append((other_key, new_item))
            eclat(prefix+[key], sorted(suffix, key=lambda item: len(item[1]), reverse=True), min_support, freq_items)
    return freq_items

def eclat_zc(data_set, min_support=1):
    # Reverse data
    data = {}
    trans_num = 0
    for trans in data_set:
        trans_num += 1
        for item in trans:
            if item not in data:
                data[item] = set()
            data[item].add(trans_num)
    print("ProcessedData",data)
    freq_items = {}
    print("ARM Start")
    freq_items = eclat([], sorted(data.items(), key=lambda item: len(item[1]), reverse=True), min_support, freq_items)
    print("ARM Finished")
    return freq_items

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

    print("healthIndex:\n", healthIndex)
    print("crcIndex:\n", crcIndex)
    print("noncrcIndex:\n", noncrcIndex)

    healthLen = len(healthIndex)
    crcLen = len(crcIndex)
    noncrcLen = len(noncrcIndex)

    print("healthLen:\n", healthLen)
    print("crcLen:\n", crcLen)
    print("noncrcLen:\n", noncrcLen)
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
                if item[column] >= 5:
                    newValue[countrow][countcolumn] = 1
                    break
            countrow += 1

    newTitle.insert(0,'index')
    value = newValue
    title = newTitle
    print("Dataset has been combined into Genus level...\n")

    # Eliminated data with all 1 or all 0, because the existence of these microorganisms is little affected by other microorganisms
    uselessColumns = []
    usefulColumns = []
    maximum = len(value[0])
    print("Start deleting useless columns...\n")
    print("Checking columns...\n")
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
    print("Useless columns are:\n", uselessColumns)

    for i in range(uselessLen-1, -1, -1):
        del title[uselessColumns[i]]
        for item in value:
            del item[uselessColumns[i]]

    print("Delete over\n")
    print("Lenth of Title:",len(title))
    print("Lenth of Value Item:", len(value[0]))
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

def rulesGenerator(frequentPatterns, minConf, rules, minLift,lenDataset):
    i = 0
    for frequentset in frequentPatterns:
        if i % 100 == 0:
            print("Rate of progress",round(i/len(frequentPatterns),2))
        if (len(frequentset) > 1):
            getRules(frequentset, frequentset, rules, frequentPatterns, minConf, minLift,lenDataset)
        i += 1

def removeStr(set, str):
    tempSet = []
    tempSet2 = []
    for elem in set:
        if (elem != str):
            tempSet.append(elem)
    tempFrozenSet = frozenset(tempSet)
    tempSet2.append(str)
    yichuSet = frozenset(tempSet2)
    return tempFrozenSet, yichuSet

def getRules(frequentset, currentset, rules, frequentPatterns, minConf, minLift,lenDataset):
    for frequentElem in currentset:
        # Frequent itemsets left after removing certain elements
        subSet, yichuSet = removeStr(currentset, frequentElem)

        try:
            a = frequentPatterns[frequentset]
        except:
            maxRecorder = 0
            for item in frequentPatterns:
                frequentsetCount = 0
                for elem in frequentset:
                    if elem in item:
                        frequentsetCount += 1
                if frequentsetCount == len(frequentset):
                    if frequentPatterns[item] > maxRecorder:
                        maxRecorder = frequentPatterns[item]
            a = maxRecorder

        try:
            b = frequentPatterns[subSet]
        except:
            maxRecorder = 0
            for item in frequentPatterns:
                subsetCount = 0
                for elem in subSet:
                    if elem in item:
                        subsetCount += 1
                if subsetCount == len(subSet):
                    if frequentPatterns[item] > maxRecorder:
                        maxRecorder = frequentPatterns[item]
            b = maxRecorder

        try:
            c = frequentPatterns[yichuSet]
        except:
            maxRecorder = 0
            for item in frequentPatterns:
                yichusetCount = 0
                for elem in yichuSet:
                    if elem in item:
                        yichusetCount += 1
                if yichusetCount == len(yichuSet):
                    if frequentPatterns[item] > maxRecorder:
                        maxRecorder = frequentPatterns[item]
            c = maxRecorder

        confidence = a / b
        lift = a / (b*(c/lenDataset))
        if ((confidence >= minConf) & (lift >= minLift)):
            flag = False
            for rule in rules:
                num = 0
                if (rule[0] == subSet and rule[1] == frequentset - subSet):
                    flag = True
                    break
                if (rule[1] == subSet and rule[0] == frequentset - subSet and rule[3] > confidence):
                    flag = True
                    break
                if (rule[1] == subSet and rule[0] == frequentset - subSet and rule[3] < confidence):
                    rules.remove(rule)
                    break
                num += 1

            if (flag == False):
                rules.append((subSet, frequentset - subSet, a/lenDataset, confidence, lift))
            if (len(subSet) >= 2):
                getRules(frequentset, subSet, rules, frequentPatterns, minConf, minLift,lenDataset)

if __name__=='__main__':
    # Dividing the datasets according to the health state of experimental samples
    healthMetadata, crcMetadata, noncrcMetadata, healthOtudata, crcOtudata, noncrcOtudata = \
            dataDevisionOtuData('crc_baxter.metadata.clean.feather', 'crc_baxter.otu_table.clean.feather')

    # Set minimum support threshold
    minSup = 0.8
    # Set minimum confidence threshold
    minConf = 0.9
    # Set minimum lift threshold
    minLift = 1

    # Preprocess the data
    dataset, title = loadDataSet(healthOtudata)

    # Random sampling 120 cases from the health dataset
    random.seed(0)
    l1 = random.sample(range(0,len(healthOtudata)),len(crcOtudata))
    print('l1',len(l1),l1)
    dataset1 = []
    for number in l1:
        dataset1.append(dataset[number])
    dataset = dataset1

    # generate frequent and candidate itemsets
    time_start1 = time.time()
    freqItems = eclat_zc(dataset, floor(minSup * len(dataset)))
    time_end1 = time.time()
    print("freqItems", len(freqItems),freqItems)

    # Delete unclosed frequent itemsets
    closedItemsDict = {}
    for item in freqItems:
        flag = True  # True: Closed Item Set, False: Unclosed Item Set
        for tempItem in freqItems:
            if tempItem == item:
                continue
            if item.issubset(tempItem):
                if freqItems[tempItem] == freqItems[item]:
                    # has a superset and the superset support is equal to the itemset support
                    flag = False
                    break
        if flag == True:
            closedItemsDict[item] = freqItems[item]
    print('number of closed frequent itemsets',len(closedItemsDict))

    # Mining positive association rules
    rules = []
    time_start2 = time.time()
    rulesGenerator(closedItemsDict, minConf, rules, minLift, len(dataset))
    time_end2 = time.time()

    # Constructing network model
    time_start3 = time.time()
    # Sort below
    # Sort rules by length
    sortedRules = []
    maxLen = 2
    for item in rules:
        if len(item[0]) + len(item[1]) > maxLen:
            maxLen = len(item[0]) + len(item[1])
    i = 2
    lsNum = 0
    names = locals()
    while i <= maxLen:
        names['a' + str(i)] = []
        lsNum += 1
        for item in rules:
            if len(item[0]) + len(item[1]) == i:
                names['a' + str(i)].append(item)
        i += 1
    # Sort the list representing each length rule by confidence level
    for i in range(2, maxLen + 1):
        names['a' + str(i)].sort(key=lambda x: x[2], reverse=True)
    # Sort items within association rules by support value
    # First invert the data
    data = {}
    trans_num = 0
    for trans in dataset:
        trans_num += 1
        for item in trans:
            if item not in data:
                data[item] = set()
            data[item].add(trans_num)
    # Sort by Support
    for i in range(2, maxLen + 1):
        for item in names['a' + str(i)]:
            for itemset in item[0:2]:
                if len(itemset) > 1:
                    ls = []
                    for number in itemset:
                        ls.append(number)
                    ls.sort(key=lambda x: data[x], reverse=True)
                    itemset = frozenset(ls)

    # save rules
    newRules = []
    for rule in rules:
        newRules.append([])
        newRules[-1].append(list(rule[0]))
        newRules[-1].append(list(rule[1]))
        newRules[-1].append(rule[2])
        newRules[-1].append(rule[3])
        newRules[-1].append(rule[4])

    uedge = []
    dedge = []
    # Reverse traversal
    for i in range(maxLen, 1, -1):
        for item in names['a' + str(i)]:
            x = list(item[0])
            y = list(item[1])
            if len(x) > 1:
                for i in range(len(x) - 1):
                    if [title[x[i]].lstrip('g__'), title[x[i + 1]].lstrip('g__')] not in uedge:
                        uedge.append([title[x[i]].lstrip('g__'), title[x[i + 1]].lstrip('g__')])
            if len(y) > 1:
                for i in range(len(y) - 1):
                    if [title[y[i]].lstrip('g__'), title[y[i + 1]].lstrip('g__')] not in uedge:
                        uedge.append([title[y[i]].lstrip('g__'), title[y[i + 1]].lstrip('g__')])
            if [title[x[0]].lstrip('g__'), title[y[0]].lstrip('g__')] not in dedge:
                dedge.append([title[x[0]].lstrip('g__'), title[y[0]].lstrip('g__')])

    print("uedge", uedge)
    print("dedge", dedge)

    # Undirected graph
    U = nx.Graph()
    # Directed graph
    D = nx.DiGraph()

    U.add_edges_from(uedge)
    D.add_edges_from(dedge)

    # Composite network
    G = nx.compose(U, D)

    # set the layout of the network
    position = nx.kamada_kawai_layout(G)

    # draw the network model
    nx.draw(G, with_labels=True, node_size=500, node_color='green', alpha=1, linewidths=3,
            font_color='black', font_weight='bold', edge_color='grey', width=1, style='dashdot',
            pos=position)
    nx.draw(D, with_labels=True, node_size=500, node_color='green', alpha=1, linewidths=3,
            font_color='black', font_weight='bold', edge_color='grey', width=1, style='solid',
            pos=position)
    time_end3 = time.time()
    plt.savefig('Baxter-MANIEA-P-NetworkModel-' + str(minSup) + '-' + str(minConf) + '-' + str(minLift) + '.png',
                dpi=500)
    plt.show()

    with open('Baxter-MANIEA-P-AssociationRules-' + str(minSup) + '-' + str(minConf) + '-' + str(minLift) + '.csv', 'w',
              newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerow(["Frequent itemset mining time", time_end1 - time_start1])
        f_csv.writerow(["Rule generation time", time_end2 - time_start2])
        f_csv.writerow(["Network construction time", time_end3 - time_start3])
        f_csv.writerow(["Source", "End", "Supp", "Conf", "Lift"])
        f_csv.writerows(newRules)

    with open('Baxter-MANIEA-P-UndirectedEdges-' + str(minSup) + '-' + str(minConf) + '-' + str(minLift) + '.csv', 'w',
              newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerows(uedge)

    with open('Baxter-MANIEA-P-DirectedEdges-' + str(minSup) + '-' + str(minConf) + '-' + str(minLift) + '.csv', 'w',
              newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerows(dedge)

    print("Frequent itemset mining time", time_end1 - time_start1)
    print("Rule generation time", time_end2 - time_start2)
    print("Network construction time", time_end3 - time_start3)

    print("Overall network indexes")
    print("Number of nodes", len(nx.degree(G)))
    print("Number of Undirected edges", U.size())
    print("Number of Positive Directed edges", D.size())

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