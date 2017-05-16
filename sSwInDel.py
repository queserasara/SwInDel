# -*- coding: utf-8 -*-
"""
Sara's script for swapping, inserting and deleting amino acid sequences in a protein. 
Written by Sara St. James to save you the agony of repetitive tasks with high error rates. Enjoy.  
sara.stjames@gmail.com
"""

import csv
import numpy as np

# start by reading the input protein : 
sFilename='input.csv' # the folder where this is run from needs to have the initial protein to be modified
sMutations='Mutations.csv' # this is a list of the modifications requested
sFilenameOut='Output.csv' # this is the output with the modified proteins
sFilename2='InputScrubbed.csv'
mData=[]
mMutations=[]
mDataNew=[]
lOutputName=[]
sOutputName=[]
lSequenceF=[]
mDeletions=[]
mInsertions=[]
# incase there are line terminators in the protein sequence




with open(sFilename, 'rU') as csvFile:  
    readData = csv.reader(csvFile, delimiter=',')
    for row in readData:
        mData.append([row[0], row[1]])
        
#if there are hard returns (\n)s this next step eliminates them

mProteinInfo=mData[0]
vProtein=mProteinInfo[1].replace('\n','')


with open(sMutations, 'rU') as csvFile: # this is the single permutation example.. upgrade for mutiple permutations later
    readData = csv.reader(csvFile, delimiter=',')
    for row in readData:
        mMutations.append([row[0]])
        mDeletions.append([row[1]])
        mInsertions.append([row[2],row[3]])

# Total Number of Mutations Requested
iNumMut=len(mMutations)-1
# read in the protein sequece to be altered
vData=mData[0]
sProteinName=vData[0]
iMutCounter=0
while iMutCounter< iNumMut:
    lOutputNameN=[]
    vLocBreaks=[]
    vLocDelBreaks=[]
    vLocDelBreaks2=[]
    vLocBreaks.append(-1)
    vLocDelBreaks.append(-1)
    vLocInsBreak=[]
    vLocInsBreak.append(-1)
    vInsSeqBreak=[]
    vInsSeqBreak.append(-1)
    lSequence=list(vProtein) # make it into a list because strings are immutable in python
    vMutation=mMutations[iMutCounter+1]
    # check for the number of mutations required
    iLenMut=len(vMutation[0])
    iFindBreak=str(vMutation[0]).find('-')
    while iFindBreak!= -1:
        vLocBreaks.append(iFindBreak)
        iLocBreak=iFindBreak
        iFindBreak=str(vMutation[0]).find('-',iLocBreak+1)
    
    vMutationSeq=vMutation[0]
    vDeletions=mDeletions[iMutCounter+1]
    vInsertions=mInsertions[iMutCounter+1]
    lOutputName.append([sProteinName,':SWAP:',vMutationSeq])
    if len(vDeletions[0])>0:
        sMutText="".join([str(x) for x in lOutputName[iMutCounter]])
        lOutputName=lOutputName[:iMutCounter]
        lOutputName.append([sMutText,':DEL:', vDeletions[0]])
    if len(vInsertions[0])>0:
        sMutText="".join([str(x) for x in lOutputName[iMutCounter]])
        lOutputName=lOutputName[:iMutCounter]
        lOutputName.append([sMutText,':INS-LOC:', vInsertions[0],':INS-SEQ:',vInsertions[1]])
        
    # determine the number of substitutions required in each required mutation
    vLocBreaks.append(iLenMut)
    iNumSwaps=len(vLocBreaks)-1    
    if iLenMut<1: 
        iNumSwaps=0
    iCounter=0
    while iCounter< iNumSwaps:
        iStart=int(vLocBreaks[iCounter])+2
        iEnd=int(vLocBreaks[iCounter+1])-1 
        iPosition=int(vMutationSeq[iStart:iEnd])-1
        iOldValue=vMutationSeq[iStart-1]
        iNewValue=vMutationSeq[iEnd]
        if lSequence[iPosition]==iOldValue:
            lSequence[iPosition]=iNewValue
        else:
            lSequence= ' E R R O R - P L E A S E C H E C K   M U T A T I O N   S E Q U E N C E  &   I N P U T '
        iCounter=iCounter+1
    # create an index for the protein to preserve original indexing
    iProteinLength=len(lSequence)
    vProteinIndex=range(1, iProteinLength+1) 
    npSequence=np.array(lSequence)
    npProteinIndex=np.array(vProteinIndex)
    # check if there are any deletions requested    
    iLenDel=len(vDeletions[0])
    if iLenDel>1:
        # determine the number of deletions        
        lDeletions=list(vDeletions[0])
        # find the - that indicate separate deletions
        iFindBreak1=str(vDeletions[0]).find('-')
        while iFindBreak1!= -1:
            vLocDelBreaks.append(iFindBreak1)
            iLocBreak=iFindBreak1
            iFindBreak1=str(vDeletions[0]).find('-',iLocBreak+1)
        vLocDelBreaks.append(iLenDel)
        
        iFindBreak2=str(vDeletions[0]).find(':')
        while iFindBreak2!= -1:
            vLocDelBreaks2.append(iFindBreak2)
            iLocBreak=iFindBreak2
            iFindBreak2=str(vDeletions[0]).find(':',iLocBreak+1)
    iNumDels=len(vLocDelBreaks2)
    iCounterDel=0
    while iCounterDel< iNumDels:
        sDeletions=vDeletions[0]
        iNum1Start=vLocDelBreaks[iCounterDel]+1
        iNum1End=vLocDelBreaks2[iCounterDel]
        iNum2Start=vLocDelBreaks2[iCounterDel]+1
        iNum2End=vLocDelBreaks[iCounterDel+1]
        iStart=int(sDeletions[iNum1Start:iNum1End])
        iEnd=int(sDeletions[iNum2Start:iNum2End])       
        vFindLow=np.where(npProteinIndex>=iStart)
        vFindHigh=np.where(npProteinIndex<=iEnd)
        vFind=np.intersect1d(vFindLow, vFindHigh)
        iStartDel=vFind[0]
        iEndDel=vFind[len(vFind)-1]
        
        # modify the protein: 
        npLow=npSequence[:iStartDel]
        npHigh=npSequence[iEndDel+1:]
        npSequence=np.concatenate((npLow,npHigh))
        # modify the protein indexing sequence
        npLowIndex=npProteinIndex[:iStartDel]       
        npHighIndex=npProteinIndex[iEndDel+1:]
        npProteinIndex=np.concatenate((npLowIndex,npHighIndex))
        
        iCounterDel=iCounterDel+1
        
# N O W   C O M P L E T E    R E Q U E S T E D    I N S E R T I O N S
    iLenIns=len(vInsertions[0])
    iLenInsSeq=len(vInsertions[1])
    if iLenIns>1:
        # determine the number of deletions        
        
        # find the - that indicate separate deletions
        iFindBreak3=str(vInsertions[0]).find('-')
        while iFindBreak3!= -1:
            vLocInsBreak.append(iFindBreak3)
            iLocBreak=iFindBreak3
            iFindBreak3=str(vInsertions[0]).find('-',iLocBreak+1)
        
        iFindBreak4=str(vInsertions[1]).find('-')
        while iFindBreak4!= -1:
            vInsSeqBreak.append(iFindBreak4)
            iLocBreak=iFindBreak4
            iFindBreak4=str(vInsertions[1]).find('-',iLocBreak+1)
        
        vLocInsBreak.append(iLenIns)
        vInsSeqBreak.append(iLenInsSeq)
    iNumIns=len(vLocInsBreak)-1
    iCounterIns=0
    while iCounterIns<iNumIns:
        sInsertions=vInsertions[0]
        sInsertSeq=vInsertions[1]
        # find location for insert- this is a single integer
        iNum1Start=vLocInsBreak[iCounterIns]+1
        iNum1End=vLocInsBreak[iCounterIns+1]
        iInsertLoc=int(sInsertions[iNum1Start:iNum1End])+1
        iSeqStart=int(vInsSeqBreak[iCounterIns])+1
        iSeqEnd=int(vInsSeqBreak[iCounterIns+1]-1)
        vInsertedSeq=sInsertSeq[iSeqStart:iSeqEnd]
        lInsertedSeq=list(vInsertedSeq)
        npInsertedSeq=np.array(lInsertedSeq)
        # find the length of the insertion value
        iLenInsert=int(vInsSeqBreak[iCounterIns+1]-vInsSeqBreak[iCounterIns])-1
        # find the location of this in the protein index (everything is referenced to positions in the original protein)
        npFindIndex1=np.where(npProteinIndex==iInsertLoc)
        iFindIndex1=int(npFindIndex1[0])        
        npDummySeq=np.zeros([iLenInsert])
        # update the Protein index sequence
        npLowIndex=npProteinIndex[:iFindIndex1]
        npHighIndex=npProteinIndex[iFindIndex1:]
        npProteinIndex0=np.concatenate((npLowIndex,npDummySeq))
        npProteinIndex=np.concatenate((npProteinIndex0,npHighIndex))
        # update the Mutated Protein Seq
        npLowSequence=npSequence[:iFindIndex1]
        npHighSequence=npSequence[iFindIndex1:]
        npSequence0=np.concatenate((npLowSequence,npInsertedSeq))
        npSequence=np.concatenate((npSequence0,npHighSequence))                
        iCounterIns=iCounterIns+1   
    lSequence=list(npSequence)
    sSequenceN="".join([str(x) for x in lSequence])     
    lSequenceF.append([sSequenceN])    
    iMutCounter=iMutCounter+1    

with open('Output.csv', 'wb') as fCSV:
    writer = csv.writer(fCSV)
    iCounter=0
    while iCounter< iNumMut:
        sOutputName="".join([str(x) for x in lOutputName[iCounter]])
        lOut0=sOutputName
        lOut1=lSequenceF[iCounter] 
        lOut2=lOut1[0]
        
        writer.writerow([lOut0,lOut2])
        iCounter=iCounter+1




