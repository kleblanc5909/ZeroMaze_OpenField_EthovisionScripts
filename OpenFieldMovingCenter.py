# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 21:28:53 2017

@author: leblanckh
"""

import pandas as pd
import numpy as np
import os
import sys
dataFolder = "/Users/leblanckh/data/KO_WT_OpenField_RawData"
resultsColumns = ["Subject", "Group", "Gender", "Timestamp", "Notes", "Time in center (%)", "Time in center while moving (%)", "Number of movements", "Average duration of movements (s)", "Total time moving (s)", "Speed while moving (cm/s)", "Average Velocity (cm/s)", "Number of movements in center"]
myDataList = []
os.chdir(dataFolder)
FileList = os.listdir(dataFolder)

def Extract_Header_Info(Filename):
    """
    Identifies the length of the header,Extracts independent variable 
    information from the header, and assigns the column names to the row that
    contains the column information in the header
    
    inputs: Filename
    Returns: NumberofHeaderRows, Subject, Genotype, Gender, Timestamp, Notes, 
    """
    NumberofHeaderRows = int(df.iloc[0,1])
    Header = df.iloc[31:NumberofHeaderRows - 3,:]
    SubjectFilters = ["Subject", "Mouse","subject", "mouse"]
    Subject =  Header[Header[0].isin(SubjectFilters)].iloc[0,1]
    GenotypeFilters = ["Genotype", "Group", "genotype", "group", "genotype/group"]
    Genotype =  Header[Header[0].isin(GenotypeFilters)].iloc[0,1]
    GenderFilters = ["Sex", "Gender", "sex", "gender"]
    Gender =  Header[Header[0].isin(GenderFilters)].iloc[0,1]
    TimestampFilters = ["Timestamp", "timestamp"]
    Timestamp =  Header[Header[0].isin(TimestampFilters)].iloc[0,1]
    NotesFilters = ["Notes", "notes"]
    Notes =  Header[Header[0].isin(NotesFilters)].iloc[0,1]
    #print ("Subject:",Subject, "Genotype:", Genotype,"Gender:", Gender,"Timestamp:", Timestamp, "Notes", Notes)
    ColumnNames = df.iloc[[NumberofHeaderRows - 2],:].values.tolist()
    df.columns = ColumnNames
    return NumberofHeaderRows,Subject,Genotype,Gender,Timestamp,Notes

for File in FileList:
    if not File.endswith('.xlsx'):
        print ("skipping file named", File)
        continue
    df = pd.read_excel(File, header = None)
    #extract information from the header in the raw data from Noldus
    
    NumberofHeaderRows,Subject,Genotype,Gender,Timestamp,Notes = Extract_Header_Info(File)
    
    #set up a data block to operate on and clean the data
    if df['Trial time'].iloc[-1] > 601:
        DataBlock = df.loc[NumberofHeaderRows +1801:,'Trial time':'Result 1']
    else:
        DataBlock = df.loc[NumberofHeaderRows +1:,'Trial time':'Result 1']
    firstRow = DataBlock[:1]
    lastRow = DataBlock[-1:]
    
    firstRow[firstRow == '-'] = 0
    lastRow[lastRow == '-'] = 0
    
    DataBlock[:1] = firstRow
    DataBlock[-1:] = lastRow
    #replace missing data (-) with NaN, then interpolate
    DataBlock.replace('-',np.NaN,inplace = True)
    DataBlock.interpolate(method = 'values', axis = 0, inplace = True)
    isMovingData = DataBlock.loc[NumberofHeaderRows +1:,'Movement(Moving / Center-point)']    
    isMovingData[isMovingData >= 0.5] = 1
    isMovingData[isMovingData < 0.5] = 0
    isInCenter = DataBlock.loc[:,'In zone']
    isInCenter[isInCenter >= 0.5] = 1
    isInCenter[isInCenter < 0.5] = 0
    Velocity = DataBlock.loc [:, 'Velocity']
    
    PercentTimeinCenter = sum (isInCenter[isInCenter == 1])/len(isInCenter)*100
    AvgVelocity = sum (Velocity)/len(Velocity)
    
    #slice the data based on the starting and ending rows for movement
    isMovingData.iat[0] = 0
    isMovingData.iat[-1] = 0
    movingTrans = isMovingData.diff()
    startingPoints = movingTrans[movingTrans == 1].index.tolist()
    endingPoints = movingTrans[movingTrans == -1].index.tolist()
    CenterTrue = 0
    inCenterCutOff = 15
    FramesCenterMoving = 0
    MovementBlocks = 0
    TotalMoveDuration = 0
    TotalVelocity = 0
    TotalFramesMoving = 0
    
    #loop over the movement data and identify if mouse enters Center
    if len(startingPoints) != len(endingPoints):
        print ("Uneven start and end pairs. There are", (len(startingPoints)), "starting points and" , (len(endingPoints)) ,"endingPoints" )
    for i in range(len(startingPoints)):
        currentStart = startingPoints[i]
        currentEnd = endingPoints[i] - 1
        #NOTE: endingPoints are one frame beyond the end of the movement.  To correct for this, I have adjusted the current end point back one frame.
        MovementBlocks += 1
        MovementDuration = (currentEnd-currentStart)/30
        TotalMoveDuration = TotalMoveDuration + MovementDuration
        FramesMoving = (currentEnd-currentStart)
        TotalFramesMoving = TotalFramesMoving + FramesMoving
        MoveVelocity = sum (Velocity.loc[currentStart:currentEnd])
        TotalVelocity = TotalVelocity + MoveVelocity
        isCenterSpan = isInCenter[currentStart:currentEnd]
        numCenter = len(isCenterSpan[isCenterSpan == 1])
        if numCenter > inCenterCutOff:
            CenterTrue +=1
            FramesCenterMoving = FramesCenterMoving + numCenter
    PercentTimeinCenterMoving = FramesCenterMoving/len(isInCenter)*100
    AvgMoveDuration = TotalMoveDuration/MovementBlocks
    AvgVelocityMoving = TotalVelocity/TotalFramesMoving

    
    
    Dataz = [Subject,Genotype, Gender, Timestamp, Notes, PercentTimeinCenter, PercentTimeinCenterMoving, MovementBlocks, AvgMoveDuration, TotalMoveDuration, AvgVelocityMoving, AvgVelocity, CenterTrue]
    myDataList.append(Dataz)
    print (myDataList)
    
resultsDF = pd.DataFrame(data = myDataList, columns=resultsColumns)
os.chdir("/Users/leblanckh/data/KO_WT_OpenField_RawData/OutputFiles")
writer = pd.ExcelWriter('KO_and_WT_OpenField_Movement_Analysis.xlsx')
resultsDF.to_excel(writer,'Sheet1')
writer.save()


