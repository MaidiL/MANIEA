import sys
import time
from numpy import *
import feather
import csv
import networkx as nx
import matplotlib.pyplot as plt
import random
type = sys.getfilesystemencoding()

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
    print(len(value[0]),len(title))
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

# Construct candidate 1-itemset C1
def createC1(dataSet):
    C1 = []
    for transaction in dataSet:
        for item in transaction:
            if not [item] in C1:
                C1.append([item])
    C1.sort()
    return list(map(frozenset, C1))

# The candidate set Ck is converted to frequent itemset Lk
# D: Original data set
# Cn: candidate set item CK
# Minsupport: threshold of support
def scanD(D, Ck, minSupport):
    # Count the number of occurrences of the candidate set
    ssCnt = {}
    for tid in D:
        for can in Ck:
            if can.issubset(tid):
                if can not in ssCnt.keys(): ssCnt[can] = 1
                else: ssCnt[can] += 1
    numItems = float(len(D))
    # Frequent itemsets LK generated by candidate itemset Cn
    Lk= []
    # Support Dictionary of candidate itemset Cn
    supportData = {}
    # The support of candidate itemsets is calculated
    for key in ssCnt:
        support = ssCnt[key] / numItems
        if support >= minSupport:
            Lk.append(key)
        supportData[key] = support
    return Lk, supportData

# The frequent Lk-1-itemsets are transformed into candidate k-itemsets by splicing
def aprioriGen(Lk_1, k):
    Ck = []
    lenLk = len(Lk_1)
    for i in range(lenLk):
        L1 = list(Lk_1[i])[:k - 2]
        L1.sort()
        for j in range(i + 1, lenLk):
            # If the first k-2 items are the same, merge the two sets
            L2 = list(Lk_1[j])[:k - 2]
            L2.sort()
            if L1 == L2:
                Ck.append(Lk_1[i] | Lk_1[j])
    return Ck

def apriori(dataSet, minSupport):
    # Construct candidate 1-itemset
    C1 = createC1(dataSet)
    # The candidate itemset Ck is converted to frequent itemset Lk
    L1, supportData = scanD(dataSet, C1, minSupport)
    L = [L1]
    k = 2
    while (len(L[k-2]) > 0):
        Lk_1 = L[k-2]
        Ck = aprioriGen(Lk_1, k)
        Lk, supK = scanD(dataSet, Ck, minSupport)
        supportData.update(supK)
        print('mining frequent',k,'-itemsets')
        L.append(Lk)
        k += 1
    return L, supportData

# Association rule generating function
# L: Frequent itemset list
# supportData: dictionary containing support data of frequent itemsets
# minConf: minimum confidence threshold
def generateRules(L, supportData, minConf=0.7):
    bigRuleList = []
    for i in range(1, len(L)):
        for freqSet in L[i]:
            H1 = [frozenset([item]) for item in freqSet]
            if (i > 1):
                rulesFromConseq(freqSet, H1, supportData, bigRuleList, minConf)
            else:
                calcConf(freqSet, H1, supportData, bigRuleList, minConf)
    return bigRuleList

# Auxiliary function -- calculate the rule's confidence conf, and filter out the rules that meet the minimum confidence requirements
# Calculate the confidence of rules and find the rules that meet the minimum confidence requirements. Function returns a list of rules that meet the minimum confidence requirement
# This list of rules is added to the bigrulelist of the main function (via the brl parameter).
# The return value prunedH holds the right part of the rule, which will be used in the next function rulesfromconseq().
def calcConf(freqSet, H, supportData, brl, minConf=0.7):
    # Evaluate the candidate rule set
    prunedH = []
    for conseq in H:
        conf = supportData[freqSet] / supportData[freqSet - conseq]
        if conf >= minConf:
            print (freqSet - conseq, '-->', conseq, 'conf:', conf)
            brl.append((freqSet - conseq, conseq, conf))
            prunedH.append(conseq)
    return prunedH

# Auxiliary function: generate the next level candidate rule set according to the current candidate rule set H
def rulesFromConseq(freqSet, H, supportData, brl, minConf=0.7):
    # Generating candidate rule sets
    m = len(H[0])
    if (len(freqSet) > (m + 1)):
        Hmpl = aprioriGen(H, m + 1)
        Hmpl = calcConf(freqSet, Hmpl, supportData, brl, minConf)
        if (len(Hmpl) > 1):
            rulesFromConseq(freqSet, Hmpl, supportData, brl, minConf)

def rulesGenerator(frequentPatterns, minConf, rules, minLift,lenDataset):
    i = 0
    for frequentset in frequentPatterns:
        if i % 10 == 0:
            print("rate of progress:",round(i/len(frequentPatterns),2))
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
    healthMetadata, crcMetadata, noncrcMetadata, healthOtudata, crcOtudata, noncrcOtudata = \
            dataDevisionOtuData('crc_baxter.metadata.clean.feather', 'crc_baxter.otu_table.clean.feather')
    # Set minimum support threshold
    minSup = 0.8
    # Set minimum confidence threshold
    minConf = 0.9
    # Set minimum lift threshold
    minLift = 1
    dataset, title = loadDataSet(healthOtudata)
    print("dataset",len(dataset),dataset)
    print("title",title)

    # Random sampling 120 cases from the health dataset
    random.seed(0)
    l1 = random.sample(range(0, len(healthOtudata)), len(crcOtudata))
    print('l1', len(l1), l1)
    dataset1 = []
    for number in l1:
        dataset1.append(dataset[number])
    dataset = dataset1
    print("dataset", len(dataset), dataset)

    timeRecorder = []
    minSupRecorder = []
    freqPatternRecorder = []
    ruleNumRecorder = []

    # mining frequent itemsets
    time_start1 = time.time()
    L, supportData = apriori(dataset, minSup)
    time_end1 = time.time()
    freqItemset = {}
    for itemset in L:
        for item in itemset:
            freqItemset[item] = round(supportData[item] * len(dataset),0)
    print("freqItemset", len(freqItemset), freqItemset)
    rules = []

    # generating association rules
    print("\nGenerating Association Rules...")
    time_start2 = time.time()
    rulesGenerator(freqItemset, minConf, rules, minLift, len(dataset))
    time_end2 = time.time()
    print("\nAssociation Rules (Min Conf:", minConf, "Min Lift:", minLift, "):")
    ruleCount = 0
    for item in rules:
        ls1 = []
        for number in item[0]:
            ls1.append(title[number].lstrip('g__'))
        ls2 = []
        for number in item[1]:
            ls2.append(title[number].lstrip('g__'))
        ruleCount += 1
    ruleNumRecorder.append(ruleCount)

    # Sort below
    # Sort rules by length
    time_start3 = time.time()
    sortedRules = []
    maxLen = 2
    for item in rules:
        if len(item[0]) + len(item[1]) > maxLen:
            maxLen = len(item[0]) + len(item[1])
    print("maxLen", maxLen)
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

    # Save Rules
    newRules= []
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
                for i in range(len(x)-1):
                    if [title[x[i]].lstrip('g__'),title[x[i+1]].lstrip('g__')] not in uedge:
                        uedge.append([title[x[i]].lstrip('g__'),title[x[i+1]].lstrip('g__')])
            if len(y) > 1:
                for i in range(len(y)-1):
                    if [title[y[i]].lstrip('g__'),title[y[i+1]].lstrip('g__')] not in uedge:
                        uedge.append([title[y[i]].lstrip('g__'),title[y[i+1]].lstrip('g__')])
            if [title[x[0]].lstrip('g__'),title[y[0]].lstrip('g__')] not in dedge:
                dedge.append([title[x[0]].lstrip('g__'),title[y[0]].lstrip('g__')])

    print("uedge",uedge)
    print("dedge",dedge)

    # Undirected graph
    U = nx.Graph()
    # Directed graph
    D = nx.DiGraph()

    U.add_edges_from(uedge)
    D.add_edges_from(dedge)

    # Composite Chart
    G = nx.compose(U,D)
    position = nx.kamada_kawai_layout(G)
    de = G.degree()
    de2 = []
    for item in de:
        de2.append(item[1] * 50)

    nx.draw(G, with_labels=True, node_size=de2, node_color='green', alpha=1, linewidths=3,
            font_color='black', font_weight='bold', edge_color='grey', width=1, style='dashdot',
            pos=position)

    nx.draw(D, with_labels=False, node_size=de2, node_color='green', alpha=1, linewidths=3,
            font_color='black', font_weight='bold', edge_color='cornflowerblue', width=1, style='solid',
            pos=position)
    time_end3 = time.time()
    plt.savefig('Baxter-Apriori-NetworkModel-' +str(minSup) +'-'+str(minConf)+ '.png', dpi=500)
    plt.show()

    with open('Baxter-Apriori-AssociationRules-'+str(minSup)+'-'+str(minConf)+'-'+str(minLift)+'.csv', 'w', newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerow(["Frequent itemset mining time", time_end1-time_start1])
        f_csv.writerow(["Rule generation time", time_end2-time_start2])
        f_csv.writerow(["Network construction time", time_end3 - time_start3])
        f_csv.writerow(["Source", "End", "Supp", "Conf", "Lift"])
        f_csv.writerows(newRules)

    with open('Baxter-Apriori-UndirectedEdges-' + str(minSup) + '-' + str(minConf) + '.csv', 'w', newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerows(uedge)

    with open('Baxter-Apriori-DirectedEdges-' + str(minSup) + '-' + str(minConf) + '.csv', 'w', newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerows(dedge)

    print("Frequent itemset mining time", time_end1 - time_start1)
    print("Rule generation time", time_end2 - time_start2)
    print("Network construction time", time_end3 - time_start3)

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