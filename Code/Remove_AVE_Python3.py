#
# Copyright 2017 Atomwise Inc.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
# associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

from scipy.spatial.distance import cdist, is_valid_dm
import argparse
import sys
from random import shuffle
import os
import numpy as np
import random
from multiprocessing import Process, Queue, Pool, Array
from scipy import stats
from time import time, asctime, localtime, strftime
from sklearn import svm
import itertools
import pickle
import gzip

from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.rdMolDescriptors import GetHashedAtomPairFingerprint, GetHashedAtomPairFingerprintAsBitVect
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.PropertyMol import PropertyMol

POP_SIZE_FACTOR = 100
NEXT_GEN_FACTOR = 20
MAX_SET_OVERLAP_FRACTION = 0.8
RM_PROB = 0.2
ADD_PROB = 0.2
MIN_AV = 10
MIN_IV = 10
MIN_AT = 30
MIN_IT = 30
FINAL_POP_OBJ = 0.01


def readMols(file):
    fileName, fileExtension = os.path.splitext(file)
    mols = []
    if(fileExtension == ".smi"):
        f = open(file, 'r')
        l = f.readline()
        f.close()
        if '\t' in l:
            mols = Chem.SmilesMolSupplier(
                file, delimiter='\t', titleLine=False)
        else:
            mols = Chem.SmilesMolSupplier(file, delimiter=' ', titleLine=False)
    elif(fileExtension == ".mol2"):
        print("mol2 file is currently not supported")
        sys.exit(1)
        mols = Chem.MolFromMol2File(file)
    elif(fileExtension == ".sdf"):
        mols = Chem.SDMolSupplier(file)
    else:
        raise Exception("un-supported input file format: " +
                        fileExtension + " . ")
    return mols


def get_fp(mols):
    fps = []
    if(args.fpType == 'ECFP4'):
        for x in mols:
            if(x):
                z = AllChem.GetMorganFingerprintAsBitVect(
                    x, 2)  # , nBits=4096 )
                fps.append(z)
    if(args.fpType == 'ECFP6'):
        for x in mols:
            if(x):
                z = AllChem.GetMorganFingerprintAsBitVect(
                    x, 3)  # , nBits=4096 )
                fps.append(z)
    if(args.fpType == 'ECFP12'):
        for x in mols:
            if(x):
                z = AllChem.GetMorganFingerprintAsBitVect(
                    x, 6)  # , nBits=4096 )
                fps.append(z)
    if(args.fpType == 'MACCS'):
        for x in mols:
            if(x):
                z = Chem.MACCSkeys.GenMACCSKeys(x)
                fps.append(z)
    if (args.fpType == 'simple'):
        describer = MUVDescriptors()
        for x in mols:
            if (x):
                z = describer.calculate_descriptors(x)
                fps.append(z)
    if(args.fpType == 'Daylight'):
        for x in mols:
            if(x):
                z = FingerprintMols.FingerprintMol(x)
                fps.append(z)
    if (args.fpType == 'AP'):
        for x in mols:
            if (x):
                z = GetHashedAtomPairFingerprintAsBitVect(x, nBits=4096)
                #z=Pairs.GetAtomPairFingerprint( x )
                fps.append(z)
    return fps


def Cdist(params):
    return cdist(params[0], params[1], params[2])


def calcDistMat(fp1, fp2, distType):
    if(args.numWorkers > 1):
        interval, reminder = divmod(len(fp1), min(len(fp1), args.numWorkers-1))
        interval = max(interval, 1)
        chuncks = pool.map(Cdist, [(fp1[r:min(
            r+interval, len(fp1))], fp2, distType) for r in range(0, len(fp1), interval)])
        return np.vstack(chuncks)
    else:
        return cdist(fp1, fp2, distType)


def checkPopSim(params):
    activeSetV, inactiveSetV, activeSetT, inactiveSetT, combinedPopJ = params[:]
    activeSetV_tmp = combinedPopJ[1][0]
    inactiveSetV_tmp = combinedPopJ[1][1]
    activeSetT_tmp = combinedPopJ[1][2]
    inactiveSetT_tmp = combinedPopJ[1][3]

    joinedActivesV = activeSetV + activeSetV_tmp
    joinedInactivesV = inactiveSetV + inactiveSetV_tmp
    joinedActivesT = activeSetT + activeSetT_tmp
    joinedInactivesT = inactiveSetT + inactiveSetT_tmp

    overlapActivesV = len(
        set([x for x in joinedActivesV if joinedActivesV.count(x) > 1]))
    overlapInactivesV = len(
        set([x for x in joinedInactivesV if joinedInactivesV.count(x) > 1]))
    overlapActivesT = len(
        set([x for x in joinedActivesT if joinedActivesT.count(x) > 1]))
    overlapInactivesT = len(
        set([x for x in joinedInactivesT if joinedInactivesT.count(x) > 1]))

    if(1.0*overlapActivesV/len(activeSetV) > MAX_SET_OVERLAP_FRACTION or
        1.0*overlapInactivesV/len(inactiveSetV) > MAX_SET_OVERLAP_FRACTION or
        1.0*overlapActivesT/len(activeSetT) > MAX_SET_OVERLAP_FRACTION or
            1.0*overlapInactivesT/len(inactiveSetT) > MAX_SET_OVERLAP_FRACTION):
        return True
    return False


def calcPopObj(params):
    cvSet, aa_D_ref, ii_D_ref, ai_D_ref = params[:]
    ia_D_ref = ai_D_ref.transpose()
    vaI, viI, taI, tiI = cvSet[:]

    # get the slices of the distance matrices
    aTest_aTrain_D = aa_D_ref[np.ix_(vaI, taI)]
    aTest_iTrain_D = ai_D_ref[np.ix_(vaI, tiI)]
    iTest_aTrain_D = ia_D_ref[np.ix_(viI, taI)]
    iTest_iTrain_D = ii_D_ref[np.ix_(viI, tiI)]

    aTest_aTrain_S = np.mean(
        [np.mean(np.any(aTest_aTrain_D < t, axis=1)) for t in np.linspace(0, 1.0, 50)])
    aTest_iTrain_S = np.mean(
        [np.mean(np.any(aTest_iTrain_D < t, axis=1)) for t in np.linspace(0, 1.0, 50)])
    iTest_iTrain_S = np.mean(
        [np.mean(np.any(iTest_iTrain_D < t, axis=1)) for t in np.linspace(0, 1.0, 50)])
    iTest_aTrain_S = np.mean(
        [np.mean(np.any(iTest_aTrain_D < t, axis=1)) for t in np.linspace(0, 1.0, 50)])

    return aTest_aTrain_S-aTest_iTrain_S + iTest_iTrain_S-iTest_aTrain_S


def writePop(pop, actives, inactives, iter):
    if(not os.path.exists(args.outDir)):
        os.makedirs(args.outDir)

    outActivesT = os.path.join(args.outDir, "actives.T.smi")
    outInactivesT = os.path.join(args.outDir, "inactives.T.smi")
    outActivesV = os.path.join(args.outDir, "actives.V.smi")
    outInactivesV = os.path.join(args.outDir, "inactives.V.smi")

    atFH = Chem.SmilesWriter(outActivesT, includeHeader=False)
    itFH = Chem.SmilesWriter(outInactivesT, includeHeader=False)
    avFH = Chem.SmilesWriter(outActivesV, includeHeader=False)
    ivFH = Chem.SmilesWriter(outInactivesV, includeHeader=False)

    for avIdx in pop[1][0]:
        avFH.write(actives[avIdx])
    for ivIdx in pop[1][1]:
        ivFH.write(inactives[ivIdx])
    for atIdx in pop[1][2]:
        atFH.write(actives[atIdx])
    for itIdx in pop[1][3]:
        itFH.write(inactives[itIdx])

    for f in os.listdir(args.outDir):
        if f.startswith("AVE_optim."):
            os.remove(os.path.join(args.outDir, f))
    open(os.path.join(args.outDir, "AVE_optim.Iter_"+str(iter)+"_obj_"+str(round(pop[0], 3)) +
                      "_size_"+str(len(pop[1][2]))+"-"+str(len(pop[1][3]))+"-"+str(len(pop[1][0]))+"-"+str(len(pop[1][1]))), "w").close()


parser = argparse.ArgumentParser(description='')
parser.add_argument('-fpType', default="ECFP4",
                    choices=['DayLight', 'ECFP4', 'ECFP6', 'ECFP12', 'AP', 'MACCS'])
parser.add_argument('-activeMols', required=True)
parser.add_argument('-inactiveMols', required=True)
parser.add_argument('-trainingToValidationRatio', default=5, type=int)
parser.add_argument('-outDir', default=".")
parser.add_argument('-numWorkers', default=1, type=int)
parser.add_argument('-statePickleFile')
parser.add_argument('-maxIter', type=int, default=300)
parser.add_argument('-maxNumMols', type=int, default=10000)
args = parser.parse_args()

startTime = time()

pool = None
checkPopSimPool = None
calcPopObjPool = None
if(args.numWorkers > 1):
    pool = Pool(processes=args.numWorkers)
    checkPopSimPool = Pool(processes=args.numWorkers)
    calcPopObjPool = Pool(processes=args.numWorkers)

pop = []
iterCount = 0
minObj = None
if(args.statePickleFile is None or not os.path.exists(args.statePickleFile)):
    actives = [m for m in readMols(args.activeMols) if m is not None]
    inactives = [m for m in readMols(args.inactiveMols) if m is not None]

    activesFP = get_fp(actives)
    inactivesFP = get_fp(inactives)

    combinedA = list(zip(actives, activesFP))
    combinedI = list(zip(inactives, inactivesFP))

    shuffle(combinedA)
    shuffle(combinedI)

    actives, activesFP = zip(*combinedA)
    inactives, inactivesFP = zip(*combinedI)

    actives = actives[:args.maxNumMols]
    inactives = inactives[:args.maxNumMols]
    activesFP = activesFP[:args.maxNumMols]
    inactivesFP = inactivesFP[:args.maxNumMols]

    print("read", len(actives), "actives and", len(inactives), "inactives")

    print("calc aa_D_ref")
    aa_D_ref = calcDistMat(activesFP, activesFP, 'jaccard')
    print("calc ii_D_ref")
    ii_D_ref = calcDistMat(inactivesFP, inactivesFP, 'jaccard')
    print("calc ai_D_ref")
    ai_D_ref = calcDistMat(activesFP, inactivesFP, 'jaccard')
    print("done")

    activeIndices = list(range(len(actives)))
    inactiveIndices = list(range(len(inactives)))

    splitActives = int(len(activeIndices)/args.trainingToValidationRatio)
    splitInactives = int(len(inactiveIndices)/args.trainingToValidationRatio)
    assert(splitActives > 0 and splitInactives > 0)
    minObj = 99999

else:
    fh = None
    if(args.statePickleFile.endswith(".gz")):
        fh = gzip.open(args.statePickleFile, "rb")
    else:
        fh = open(args.statePickleFile, "rb")
    actives, inactives, aa_D_ref, ii_D_ref, ai_D_ref, iterCount, pop, minObj = pickle.load(
        fh)

# randomly select a population of CV sets
POP_SIZE = POP_SIZE_FACTOR
while(len(pop) < POP_SIZE):
    shuffle(activeIndices)
    shuffle(inactiveIndices)
    pop.append((activeIndices[:splitActives], inactiveIndices[:splitInactives],
                activeIndices[splitActives:], inactiveIndices[splitInactives:]))

while(iterCount < args.maxIter):
    iterCount += 1

    print("calculate objectives for the population")
    objs = []
    if(args.numWorkers > 1):
        objs = calcPopObjPool.map(
            calcPopObj, [(cvSet, aa_D_ref, ii_D_ref, ai_D_ref) for cvSet in pop])
    else:
        for cvSet in pop:
            objs.append(calcPopObj([cvSet, aa_D_ref, ii_D_ref, ai_D_ref]))

    # negeate the values if we are maximizing
    combinedPop = sorted(zip(objs, pop))
    if(combinedPop[0][0] < 0):
        objs = [-i for i in objs]
        combinedPop = sorted(zip(objs, pop))

    # enforce some diversity so the CV sets would not all be too similar
    print("remove similar sets")
    skipIndices = {}
    for i in range(len(combinedPop)):
        if(i in skipIndices):
            continue
        activeSetV = combinedPop[i][1][0]
        inactiveSetV = combinedPop[i][1][1]
        activeSetT = combinedPop[i][1][2]
        inactiveSetT = combinedPop[i][1][3]

        if(args.numWorkers > 1):
            results = checkPopSimPool.map(checkPopSim,
                                          [(activeSetV, inactiveSetV, activeSetT, inactiveSetT, combinedPop[j])
                                           for j in range(i+1, len(combinedPop)) if j not in skipIndices])
            for rIdx, val in enumerate(results):
                if(val):
                    skipIndices[rIdx] = 1
                else:
                    for j in range(i+1, len(combinedPop)):
                        if(j in skipIndices):
                            continue
                        val = checkPopSim(
                            (activeSetV, inactiveSetV, activeSetT, inactiveSetT, combinedPop[j]))
                        if(val):
                            skipIndices[rIdx] = 1

        numSkipped = len([x for x in skipIndices if x < i])
        if(i - numSkipped > NEXT_GEN_FACTOR):
            break

    # now, remove the similar entries
    print("removing", min(len(skipIndices.keys()),
                          NEXT_GEN_FACTOR), "similar sets")
    for i in reversed(np.sort(list(skipIndices.keys()))):
        del combinedPop[i]
        if(len(combinedPop) <= NEXT_GEN_FACTOR):
            break

    # select the top K sets
    print("population size after similarity filter: ", len(combinedPop))
    print("select the next generation")
    newPop = combinedPop[:NEXT_GEN_FACTOR]
    finalPop = combinedPop[0]

    fullPopObj = np.mean([x[0] for x in combinedPop])
    topPopObj = np.mean([x[0] for x in newPop])
    finalPopObj = finalPop[0]

    print("iter=", iterCount, "fullPopObj=", round(fullPopObj, 3), "topPopObj=", round(
        topPopObj, 3), "finalPopObj=", round(finalPopObj, 3), "minObj=", round(minObj, 3))

    if(abs(finalPopObj) < abs(minObj)):
        minObj = finalPopObj

        if(args.statePickleFile != None):
            if(not os.path.isdir(os.path.dirname(args.statePickleFile))):
                os.makedirs(os.path.dirname(args.statePickleFile))
            if(args.statePickleFile.endswith(".gz")):
                pickle.dump(([PropertyMol(m) for m in actives], [PropertyMol(m) for m in inactives], aa_D_ref, ii_D_ref,
                             ai_D_ref, iterCount, pop, minObj), gzip.open(args.statePickleFile, "wb"), protocol=pickle.HIGHEST_PROTOCOL)
            else:
                pickle.dump(([PropertyMol(m) for m in actives], [PropertyMol(m) for m in inactives], aa_D_ref, ii_D_ref,
                             ai_D_ref, iterCount, pop, minObj), open(args.statePickleFile, "wb"), protocol=pickle.HIGHEST_PROTOCOL)
        writePop(finalPop, actives, inactives, iterCount)

    if(abs(finalPopObj) < FINAL_POP_OBJ):
        break

    # breed
    print("breed")
    pop = []
    while(len(pop) < POP_SIZE):
        # randomly choose a pair
        pair = random.sample(newPop, 2)

        # for each subset of indices, select a new set from the union of indices for that subset
        newActiveIndicesV = pair[0][1][0] + pair[1][1][0]
        newInactiveIndicesV = pair[0][1][1] + pair[1][1][1]
        newActiveIndicesT = pair[0][1][2] + pair[1][1][2]
        newInactiveIndicesT = pair[0][1][3] + pair[1][1][3]

        avSize = int(len(newActiveIndicesV)/2)
        ivSize = int(len(newInactiveIndicesV)/2)
        atSize = int(len(newActiveIndicesT)/2)
        itSize = int(len(newInactiveIndicesT)/2)

        newActiveIndicesV = np.unique(newActiveIndicesV)
        newInactiveIndicesV = np.unique(newInactiveIndicesV)
        newActiveIndicesT = np.unique(newActiveIndicesT)
        newInactiveIndicesT = np.unique(newInactiveIndicesT)

        # make sure there are no overlapping molecules between training/validation sets
        newActiveIndices = list(
            np.hstack((newActiveIndicesV, newActiveIndicesT)))
        newInactiveIndices = list(
            np.hstack((newInactiveIndicesV, newInactiveIndicesT)))
        overlapActives = list(
            set([x for x in newActiveIndices if newActiveIndices.count(x) > 1]))
        overlapInactives = list(
            set([x for x in newInactiveIndices if newInactiveIndices.count(x) > 1]))
        shuffle(overlapActives)
        shuffle(overlapInactives)
        for idx, overlapA in enumerate(overlapActives):
            if(idx % 2):
                newActiveIndicesV = np.delete(
                    newActiveIndicesV, np.where(newActiveIndicesV == overlapA))
            else:
                newActiveIndicesT = np.delete(
                    newActiveIndicesT, np.where(newActiveIndicesT == overlapA))
        for idx, overlapI in enumerate(overlapInactives):
            if(idx % 2):
                newInactiveIndicesV = np.delete(
                    newInactiveIndicesV, np.where(newInactiveIndicesV == overlapI))
            else:
                newInactiveIndicesT = np.delete(
                    newInactiveIndicesT, np.where(newInactiveIndicesT == overlapI))

        newActiveIndices = list(
            np.hstack((newActiveIndicesV, newActiveIndicesT)))
        newInactiveIndices = list(
            np.hstack((newInactiveIndicesV, newInactiveIndicesT)))
        assert(
            len(set([x for x in newActiveIndices if newActiveIndices.count(x) > 1])) == 0)
        assert(
            len(set([x for x in newInactiveIndices if newInactiveIndices.count(x) > 1])) == 0)

        avSize = min(avSize, len(newActiveIndicesV))
        ivSize = min(ivSize, len(newInactiveIndicesV))
        atSize = min(atSize, len(newActiveIndicesT))
        itSize = min(itSize, len(newInactiveIndicesT))

        if(np.random.rand() < RM_PROB and avSize > MIN_AV):
            avSize -= 1
        if(np.random.rand() < ADD_PROB):
            avSize += 1
        if(np.random.rand() < RM_PROB and ivSize > MIN_IV):
            ivSize -= 1
        if(np.random.rand() < ADD_PROB):
            ivSize += 1
        if(np.random.rand() < RM_PROB and atSize > MIN_AT):
            atSize -= 1
        if(np.random.rand() < ADD_PROB):
            atSize += 1
        if(np.random.rand() < RM_PROB and itSize > MIN_IT):
            itSize -= 1
        if(np.random.rand() < ADD_PROB):
            itSize += 1

        avSamp = random.sample(list(newActiveIndicesV), min(
            len(newActiveIndicesV), max(avSize, MIN_AV)))
        ivSamp = random.sample(list(newInactiveIndicesV), min(
            len(newInactiveIndicesV), max(ivSize, MIN_IV)))
        atSamp = random.sample(list(newActiveIndicesT), min(
            len(newActiveIndicesT), max(atSize, MIN_AT)))
        itSamp = random.sample(list(newInactiveIndicesT), min(
            len(newInactiveIndicesT), max(itSize, MIN_IT)))

        pop.append((avSamp, ivSamp, atSamp, itSamp))


print("Done. Total running time: %.3f sec\n" % (time() - startTime))
