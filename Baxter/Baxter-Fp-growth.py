import time
import sys
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

# Convert the list to frozenset
def transfer2FrozenDataSet(dataSet):
    frozenDataSet = {}
    for elem in dataSet:
        frozenDataSet[frozenset(elem)] = 1
    return frozenDataSet

# Classes of Fp-Tree node
class TreeNode:
    def __init__(self, nodeName, count, nodeParent):
        self.nodeName = nodeName
        self.count = count
        self.nodeParent = nodeParent
        self.nextSimilarItem = None
        self.children = {}

    def increaseC(self, count):
        self.count += count

    def disp(self, ind=1):
        print ('  '*ind, self.nodeName, ' ', self.count)
        for child in self.children.values():
            child.disp(ind+1)

def createFPTree(frozenDataSet, minSupport):
    # First let's think about the process of creating an FP tree
    '''
    Usually we have two steps:
    First step: Scan the dataset once, filter out the infrequent items, get a frequent set of items, and store the number of times they occur, because the second time you scan a record, you need to
    Step 2: Scan the dataset again, remove the infrequent items for each record, sort the items in the record from the highest to the lowest, and insert the FP tree
    '''
    # Step 1
    headPointTable = {}
    for items in frozenDataSet:
        for item in items:
            headPointTable[item] = headPointTable.get(item, 0) + frozenDataSet[items]
    headPointTable = {k:v for k,v in headPointTable.items() if v >= minSupport}
    frequentItems = set(headPointTable.keys())
    if len(frequentItems) == 0: return None, None

    for k in headPointTable:
        headPointTable[k] = [headPointTable[k], None]
    fptree = TreeNode("null", 1, None)
    #scan dataset at the second time, filter out items for each record
    for items,count in frozenDataSet.items():
        frequentItemsInRecord = {}
        for item in items:
            if item in frequentItems:
                frequentItemsInRecord[item] = headPointTable[item][0]
        if len(frequentItemsInRecord) > 0:
            orderedFrequentItems = [v[0] for v in sorted(frequentItemsInRecord.items(), key=lambda v:v[1], reverse = True)]
            updateFPTree(fptree, orderedFrequentItems, headPointTable, count)
    return fptree, headPointTable

def updateFPTree(fptree, orderedFrequentItems, headPointTable, count):
    #handle the first item
    if orderedFrequentItems[0] in fptree.children:
        fptree.children[orderedFrequentItems[0]].increaseC(count)
    else:
        fptree.children[orderedFrequentItems[0]] = TreeNode(orderedFrequentItems[0], count, fptree)

        #update headPointTable
        if headPointTable[orderedFrequentItems[0]][1] == None:
            headPointTable[orderedFrequentItems[0]][1] = fptree.children[orderedFrequentItems[0]]
        else:
            updateHeadPointTable(headPointTable[orderedFrequentItems[0]][1], fptree.children[orderedFrequentItems[0]])
    #handle other items except the first item
    if(len(orderedFrequentItems) > 1):
        updateFPTree(fptree.children[orderedFrequentItems[0]], orderedFrequentItems[1::], headPointTable, count)

def updateHeadPointTable(headPointBeginNode, targetNode):
    while(headPointBeginNode.nextSimilarItem != None):
        headPointBeginNode = headPointBeginNode.nextSimilarItem
    headPointBeginNode.nextSimilarItem = targetNode

def mineFPTree(headPointTable, prefix, frequentPatterns, minSupport):
    #for each item in headPointTable, find conditional prefix path, create conditional fptree, then iterate until there is only one element in conditional fptree
    headPointItems = [v[0] for v in sorted(headPointTable.items(), key = lambda v:v[1][0])]
    if(len(headPointItems) == 0): return

    for headPointItem in headPointItems:
        newPrefix = prefix.copy()
        newPrefix.add(headPointItem)
        support = headPointTable[headPointItem][0]
        frequentPatterns[frozenset(newPrefix)] = support

        prefixPath = getPrefixPath(headPointTable, headPointItem)
        if(prefixPath != {}):
            conditionalFPtree, conditionalHeadPointTable = createFPTree(prefixPath, minSupport)
            if conditionalHeadPointTable != None:
                mineFPTree(conditionalHeadPointTable, newPrefix, frequentPatterns, minSupport)

def getPrefixPath(headPointTable, headPointItem):
    prefixPath = {}
    beginNode = headPointTable[headPointItem][1]
    prefixs = ascendTree(beginNode)
    if((prefixs != [])):
        prefixPath[frozenset(prefixs)] = beginNode.count
    while(beginNode.nextSimilarItem != None):
        beginNode = beginNode.nextSimilarItem
        prefixs = ascendTree(beginNode)
        if (prefixs != []):
            prefixPath[frozenset(prefixs)] = beginNode.count
    return prefixPath

# Find the path from the current node to the root node (bottom-up), which is stored in prefixPath with the name of the item
def ascendTree(treeNode):
    prefixs = []
    while((treeNode.nodeParent != None) and (treeNode.nodeParent.nodeName != 'null')):
        treeNode = treeNode.nodeParent
        prefixs.append(treeNode.nodeName)
    return prefixs

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

    timeRecorder = []
    minSupRecorder = []
    freqPatternRecorder = []
    ruleNumRecorder = []

    # Calculate minimum occureence frequency
    minSupport = floor(minSup * len(dataset))
    frozenDataSet = transfer2FrozenDataSet(dataset)
    frequentPatterns = {}
    prefix = set([])
    time_start1 = time.time()
    fptree, headPointTable = createFPTree(frozenDataSet, minSupport)
    mineFPTree(headPointTable, prefix, frequentPatterns, minSupport)
    time_end1 = time.time()

    rules = []
    print('frequentPatterns',len(frequentPatterns),frequentPatterns)
    dic1SortList = sorted(frequentPatterns.items(), key=lambda x: x[1], reverse=True)
    print("\nGenerating Association Rules...")
    time_start2 = time.time()
    rulesGenerator(frequentPatterns, minConf, rules, minLift, len(dataset))
    time_end2 = time.time()
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

    # #Sort the list representing each length rule by confidence level
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
    nx.draw(D, with_labels=True, node_size=de2, node_color='green', alpha=1, linewidths=3,
            font_color='black', font_weight='bold', edge_color='cornflowerblue', width=1, style='solid',
            pos=position)
    time_end3 = time.time()
    plt.savefig('Baxter-Fpgrowth-NetworkModel-'+str(minSup)+'-'+str(minConf)+'-'+str(minLift)+'.png',dpi=500)
    plt.show()

    with open('Baxter-Fpgrowth-AssociationRules-'+str(minSup)+'-'+str(minConf)+'-'+str(minLift)+'.csv', 'w', newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerow(["Frequent itemset mining time", time_end1-time_start1])
        f_csv.writerow(["Rule generation time", time_end2-time_start2])
        f_csv.writerow(["Network construction time", time_end3 - time_start3])
        f_csv.writerow(["Source", "End", "Supp", "Conf", "Lift"])
        f_csv.writerows(newRules)

    with open('Baxter-Fpgrowth-UndirectedEdges-' + str(minSup) + '-' + str(minConf) + '-' + str(minLift) + '.csv', 'w', newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerows(uedge)

    with open('Baxter-Fpgrowth-DirectedEdges-' + str(minSup) + '-' + str(minConf) + '-' + str(minLift)  + '.csv', 'w', newline='') as f:
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