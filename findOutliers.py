""" This code finds outliers using the outlier scores
    by computing mean and variance of the outlier scores """

import collections

inputFile = "outlierScores"
outputFile = "outliers"

motifType = {}
def readMotifTypes():
    global motifType
    fin = open("reSampledMotifs", "r")
    for i,line in enumerate(fin):
        line = line.rstrip().split()
        motifType[i] = int(line[4])
    fin.close()

delta = 1.0
def findOutliers():
    """ Reads outlier scores, computes mean and variance and then finds outliers """

    muDict = collections.defaultdict(float)
    countDict = collections.defaultdict(int)

    triples = []
    totalOutlierScore = 0.0
    fin = open(inputFile, "r")
    for line in fin:
        triple = line.rstrip().split()
        mType = motifType[int(triple[0])]
        muDict[mType]+=float(triple[2])
        countDict[mType]+=1
        triples.append(triple)
    fin.close()

    for mType in muDict:
        muDict[mType]/=countDict[mType]


    sigmaDict = collections.defaultdict(float)
    sigma = 0.0
    for triple in triples:
        mType = motifType[int(triple[0])]
        x = float(triple[2])
        sigmaDict[mType] += (x-muDict[mType])**2

    for mType in sigmaDict:
        sigmaDict[mType]/=countDict[mType]
        sigmaDict[mType] = sigmaDict[mType]**(0.5)

    print len(triples)

    ## 23 : [10, 100] means 23 nbrs wale 100 mein se 10 are outliers
    neighbourOutlierDict = {}
    for mType in muDict:
        neighbourOutlierDict[mType] = {}

    fout = open(outputFile, "w")
    for triple in triples:
        mType = motifType[int(triple[0])]
        if int(triple[1]) not in neighbourOutlierDict[mType]:
            neighbourOutlierDict[mType][int(triple[1])] = [0, 0]

        neighbourOutlierDict[mType][int(triple[1])][1] += 1
        x = float(triple[2])
        ## An outlier
        threshold = muDict[mType] + delta * sigmaDict[mType]
        if x > threshold:
            fout.write(str(mType) + " " + triple[1] + " " +triple[2] + "\n")
            neighbourOutlierDict[mType][int(triple[1])][0] += 1
    fout.close()


    fout = open("results.csv", "w")
    for mType in muDict:
        print mType
        fout.write(str(mType) + "\n")
        keylist = neighbourOutlierDict[mType].keys()
        # keylist.sort()
        for key in range(1,51):
            nbrCount = key
            if nbrCount not in keylist:
                print nbrCount, 0, 0, 0
                fout.write(str(nbrCount) + " " + str(0) + " " + str(0) + " " + str(0.0) + "\n")
                continue

            outlierCount = neighbourOutlierDict[mType][key][0]
            totalMotifs = neighbourOutlierDict[mType][key][1]
            # print nbrCount, float(outlierCount)/float(totalMotifs)
            print nbrCount, outlierCount, totalMotifs, float(outlierCount)/float(totalMotifs)
            fout.write(str(nbrCount) + " " + str(outlierCount) + " " + str(totalMotifs) + " " + str(float(outlierCount)/float(totalMotifs)) + "\n")
    fout.close()


if __name__ == "__main__":
    readMotifTypes()
    findOutliers()
