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

# Candidate Itemset Global Variable to hold the set of candidate itemsets generated during the mining of frequent itemsets
candidateRecorder = {}

def eclat2(prefix, items, min_support, freq_items):
    global candidateRecorder
    while items:
        # get itemset and corresponding transaction number from items one by one
        key, item = items.pop()
        # The number of transactions is the support of the itemset
        key_support = len(item)
        # Record Candidate Itemset
        candidateRecorder[frozenset(sorted(prefix + [key]))] = key_support
        # Determine if itemsets are frequent
        if key_support >= min_support:
            # If the itemset is frequent, add the itemset and its support to the frequent itemset set freq_items
            freq_items[frozenset(sorted(prefix + [key]))] = key_support
            suffix = []
            for other_key, other_item in items:
                new_item = item & other_item
                # Record Candidate Itemset
                candidateRecorder[frozenset(sorted(prefix + [key] + [other_key]))] = len(new_item)
                # Determine if the candidate itemset occur frequently
                # If frequent, add it to suffix
                if len(new_item) >= min_support:
                    suffix.append((other_key, new_item))
            eclat2(prefix + [key], sorted(suffix, key=lambda item: len(item[1]), reverse=True), min_support, freq_items)
    return freq_items

def eclat_zc2(data_set, min_support=1):
    # Reverse data
    data = {}
    trans_num = 0
    for trans in data_set:
        trans_num += 1
        for item in trans:
            if item not in data:
                data[item] = set()
            data[item].add(trans_num)
    print("ProcessedData", data)
    freq_items = {}
    print("ARM Start")
    freq_items = eclat2([], sorted(data.items(), key=lambda item: len(item[1]), reverse=True), min_support, freq_items)
    print("ARM Finished")
    return freq_items

# Divide the data by health status of samples
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

    print("healthIndex:\n", healthIndex)
    print("eddIndex:\n", eddIndex)

    healthLen = len(healthIndex)
    eddLen = len(eddIndex)

    print("healthLen:\n", healthLen)
    print("eddLen:\n", eddLen)
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

    print('eddDropIndex',eddDropIndex)
    print('healthDropIndex',healthDropIndex)

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

    newTitle.insert(0, 'index')
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
    for j in range(1, maximum):
        if j not in usefulColumns:
            uselessColumns.append(j)
    uselessColumns.sort()
    uselessLen = len(uselessColumns)
    print("Useless columns are:\n", uselessColumns)

    for i in range(uselessLen - 1, -1, -1):
        del title[uselessColumns[i]]
        for item in value:
            del item[uselessColumns[i]]

    print("Delete over\n")
    print("Lenth of Title:", len(title))
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
    return data, title

def rulesGenerator(frequentPatterns, minConf, rules, minLift, lenDataset):
    i = 0
    for frequentset in frequentPatterns:
        if i % 100 == 0:
            print("进度：", round(i / len(frequentPatterns), 2))
        if (len(frequentset) > 1):
            getRules(frequentset, frequentset, rules, frequentPatterns, minConf, minLift, lenDataset)
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

def getRules(frequentset, currentset, rules, frequentPatterns, minConf, minLift, lenDataset):
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
        lift = a / (b * (c / lenDataset))
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
                rules.append((subSet, frequentset - subSet, a / lenDataset, confidence, lift))
            if (len(subSet) >= 2):
                getRules(frequentset, subSet, rules, frequentPatterns, minConf, minLift, lenDataset)

if __name__ == '__main__':
    # Dividing the datasets according to the health state of experimental samples
    healthMetadata, eddMetadata, healthOtudata, eddOtudata = \
        dataDevisionOtuData('edd_singh.metadata.clean.feather', 'edd_singh.otu_table.clean.feather')

    # Set minimum support threshold
    minSup = 0.65
    # Set minimum confidence threshold
    minConf = 0.75
    # Set minimum pearson correlation coefficient threshold
    minCor = 0.1
    # Set minimum lift threshold
    minLift = 1

    # Preprocess the data
    dataset, title = loadDataSet(healthOtudata)
    print('dataset', dataset)
    print("title", title)


    # generate frequent and candidate itemsets
    time_start1 = time.time()
    freqItems = eclat_zc2(dataset, floor(minSup * len(dataset)))
    time_end1 = time.time()
    print("freqItems", len(freqItems),freqItems)
    print("candidateRecorder", len(candidateRecorder), candidateRecorder)

    # Delete unclosed candidate itemsets
    closedItemsDict = {}
    for item in candidateRecorder:
        flag = True  # True: Closed Item Set, False: Unclosed Item Set
        for tempItem in candidateRecorder:
            if tempItem == item:
                continue
            if item.issubset(tempItem):
                if candidateRecorder[tempItem] == candidateRecorder[item]:
                    # has a superset and the superset support is equal to the itemset support
                    flag = False
                    break
        if flag == True:
            closedItemsDict[item] = candidateRecorder[item]

    print('closedItemsDict',len(closedItemsDict),closedItemsDict)

    # Mining positive and negative association rules
    iterCount = 0
    dataLen = len(dataset)
    print("Mining positive and negative association rules")
    posRules = []
    negRules = []
    time_start2 = time.time()
    for itemset in closedItemsDict.keys():
        iterCount += 1
        if iterCount % 50 == 0:
            print("Association Rule Mining Process：", round(100 * iterCount / len(closedItemsDict.keys()), 1), "%")
            time_end = time.time()
            print("Time consuming：", round((time_end - time_start2) // 60, 0), "min", round((time_end - time_start2) % 60, 2),
                  "s")
        # Traverse through all possible combinations of items in candidate itemsets
        for i in range(1, len(itemset)):
            for j in combinations(itemset, i):
                # ls1: antecedent
                # ls2: consequent
                # Potential association rule : ls1 -> ls2, ls1 -> (neg) ls2, (neg) ls1 -> ls2, (neg) ls1 -> (neg) ls2
                ls1 = list(j)
                lsLeft = list(set(itemset) - (set(j)))
                ls2 = list(lsLeft)
                # calculate the pearson correlation coefficient between the occurence vectors of antecedent and consequent
                # the occurence vectors represent whether the itemsets appear in each data record
                x = []
                y = []
                for data in dataset:
                    if set(ls1) < set(data):
                        x.append(1)
                    else:
                        x.append(0)
                    if set(ls2) < set(data):
                        y.append(1)
                    else:
                        y.append(0)
                x = pd.Series(x)
                y = pd.Series(y)
                corValue = x.corr(y, method='pearson')

                if corValue >= minCor:
                    nxnyCount = 0
                    xnyCount = 0
                    nxyCount = 0
                    xyCount = 0
                    for i in range(dataLen):
                        if x[i] == 0 and y[i] == 0:
                            nxnyCount += 1
                        elif x[i] == 0 and y[i] == 1:
                            nxyCount += 1
                        elif x[i] == 1 and y[i] == 0:
                            xnyCount += 1
                        elif x[i] == 1 and y[i] == 1:
                            xyCount += 1
                    xCount = sum(x)
                    yCount = sum(y)
                    nxCount = dataLen - xCount
                    nyCount = dataLen - yCount
                    # analyze the interst indicators of the rule x->y
                    if closedItemsDict[itemset] / dataLen >= minSup:
                        if closedItemsDict[itemset] / xCount >= minConf and closedItemsDict[itemset] / xCount / (
                                yCount / dataLen) >= minLift:
                            posRules.append([ls1, ls2, round(corValue, 3), round(xyCount / dataLen, 3),
                                             round(xyCount / xCount, 3),
                                             closedItemsDict[itemset] / xCount / (yCount / dataLen), 'x->y'])
                    # analyze the interst indicators of the rule nx->ny
                    elif nxnyCount / nxCount >= minConf and nxnyCount / nxCount / (nyCount / dataLen) >= minLift:
                        posRules.append([ls1, ls2, round(corValue, 3), round(nxnyCount / dataLen, 3),
                                         round(nxnyCount / nxCount, 3), nxnyCount / nxCount / (nyCount / dataLen),
                                         'nx->ny'])

                elif corValue <= -1 * minCor:
                    nxnyCount = 0
                    xnyCount = 0
                    nxyCount = 0
                    xyCount = 0
                    for i in range(dataLen):
                        if x[i] == 0 and y[i] == 0:
                            nxnyCount += 1
                        elif x[i] == 0 and y[i] == 1:
                            nxyCount += 1
                        elif x[i] == 1 and y[i] == 0:
                            xnyCount += 1
                        elif x[i] == 1 and y[i] == 1:
                            xyCount += 1

                    xCount = sum(x)
                    yCount = sum(y)
                    nxCount = dataLen - xCount
                    nyCount = dataLen - yCount
                    # analyze the interst indicators of the rule x->ny
                    if xnyCount / xCount >= minConf and xnyCount / xCount / (nyCount / dataLen) >= minLift:
                        negRules.append([ls1, ls2, round(corValue, 3), round(xnyCount / dataLen, 3),
                                             round(xnyCount / xCount, 3), xnyCount / xCount / (nyCount / dataLen),
                                             'x->ny'])
                    # analyze the interst indicators of the rule nx->y
                    if nxyCount / nxCount >= minConf and nxyCount / nxCount / (yCount / dataLen) >= minLift:
                        negRules.append([ls1, ls2, round(corValue, 3), round(nxyCount / dataLen, 3),
                                             round(nxyCount / nxCount, 3), nxyCount / nxCount / (yCount / dataLen),
                                             'nx->y'])
    time_end2 = time.time()

    posRules2 = []
    negRules2 = []

    # Remove duplicates in positive rules and negative rules
    for item in posRules:
        if item not in posRules2:
            posRules2.append(item)
    posRules = posRules2
    for item in negRules:
        if item not in negRules2:
            negRules2.append(item)
    negRules = negRules2

    print("posRules", len(posRules), posRules)
    print("negRules", len(negRules), negRules)


    # Constructing network model
    time_start3 = time.time()
    uedge = []
    posDedge = []
    negDedge = []
    # Processing positive and negative rules separately
    # Sort below
    # Sort rules by length
    for rulePosNeg in range(0, 2):
        if rulePosNeg == 0:
            print("Processing positive rules")
            rules = posRules
        elif rulePosNeg == 1:
            print("Processing negative rules")
            rules = negRules
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
        if rulePosNeg == 0:
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
                    if [title[x[0]].lstrip('g__'), title[y[0]].lstrip('g__')] not in posDedge:
                        posDedge.append([title[x[0]].lstrip('g__'), title[y[0]].lstrip('g__')])
        if rulePosNeg == 1:
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
                    if [title[x[0]].lstrip('g__'), title[y[0]].lstrip('g__')] not in negDedge:
                        negDedge.append([title[x[0]].lstrip('g__'), title[y[0]].lstrip('g__')])

    # Undirected graph
    U = nx.Graph()
    # Posotive Directed graph
    Dpos = nx.DiGraph()
    # Negative Directed graph
    Dneg = nx.DiGraph()

    U.add_edges_from(uedge)
    Dpos.add_edges_from(posDedge)
    Dneg.add_edges_from(negDedge)

    print("Number of undirected edges", len(U.edges))
    print("Number of directed edges", len(Dpos.edges) + len(Dneg.edges))

    # Composite network
    D = nx.compose(Dpos, Dneg)
    G = nx.compose(U, D)

    # set the size of the nodes
    de = G.degree()
    de2 = []
    for item in de:
        de2.append(item[1] * 50)
    de3 = D.degree()
    de4 = []
    for item in de3:
        de4.append(item[1] * 50)

    print("G.nodes", G.nodes())
    print("G.degree", G.degree())
    print("G.edges", len(G.edges), G.edges())
    print("Dpos.nodes", Dpos.nodes())
    print("Dpos.degree", Dpos.degree())
    print("Dpos.edges", len(Dpos.edges), Dpos.edges())
    print("Dneg.nodes", Dneg.nodes())
    print("Dneg.degree", Dneg.degree())
    print("Dneg.edges", len(Dneg.edges), Dneg.edges())

    # set the layout of the network
    position = nx.kamada_kawai_layout(G)
    print("position", position)

    # draw the network model
    nx.draw(G, with_labels=True, node_size=de2, node_color='green', alpha=1, linewidths=3,
            font_color='black', edge_color='grey', width=1, style='dashdot',
            pos=position, font_size=10)
    nx.draw(Dpos, with_labels=False, node_size=de4, node_color='green', alpha=1, linewidths=3,
            font_color='black', edge_color='grey', width=1, style='solid',
            pos=position)
    nx.draw(Dneg, with_labels=False, node_size=de4, node_color='green', alpha=1, linewidths=3,
            font_color='black', edge_color='red', width=1.5, style='dashed',
            pos=position)
    time_end3 = time.time()

    plt.savefig('Singh-MANIEA-NetworkModel-' + str(minSup) + '-' + str(minConf) + '-' + str(minCor) + '.png',
                dpi=500)
    plt.show()

    with open('Singh-MANIEA-AssociationRules-' + str(minSup) + '-' + str(minConf) + '-' + str(minCor) + '.csv', 'w',
              newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerow(["Frequent itemset mining time", time_end1 - time_start1])
        f_csv.writerow(["Rule generation time", time_end2 - time_start2])
        f_csv.writerow(["Network construction time", time_end3 - time_start3])
        f_csv.writerow(["Source", "End", "Cor", "Supp", "Conf", "Lift", "Type"])
        f_csv.writerows(negRules)
        f_csv.writerows(posRules)

    with open('Singh-MANIEA-UndirectedEdges-' + str(minSup) + '-' + str(minConf) + '-' + str(minCor) + '.csv', 'w',
              newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerows(uedge)

    with open('Singh-MANIEA-PositiveDirectedEdges-' + str(minSup) + '-' + str(minConf) + '-' + str(minCor) + '.csv',
              'w', newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerows(posDedge)

    with open('Singh-MANIEA-NegativeDirectedEdges-' + str(minSup) + '-' + str(minConf) + '-' + str(minCor) + '.csv',
              'w', newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerows(negDedge)

    print("Frequent itemset mining time", time_end1 - time_start1)
    print("Rule generation time", time_end2 - time_start2)
    print("Network construction time", time_end3 - time_start3)

    print("0. Overall network indexes")
    print("Number of nodes", len(nx.degree(G)))
    print("Number of Undirected edges", U.size())
    print("Number of positive directed edges", Dpos.size())
    print("Number of negative directed edges", Dneg.size())

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