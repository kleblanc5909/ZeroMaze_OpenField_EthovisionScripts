# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 21:28:53 2017

@author: leblanckh
"""

import pandas as pd
import numpy as np
import os
import sys
dataFolder = "/Users/leblanckh/data/DREADD_ZeroMaze_RawData"
resultsColumns = ["Subject", "Drug", "Diet", "Timestamp", "Notes", "Time in open (%)", "Time in open while moving (%)", "Number of movements", "Average duration of movements (s)", "Total time moving (s)", "Speed while moving (cm/s)", "Average velocity (cm/s)", "Number of movements in open"]
myDataList = []
os.chdir(dataFolder)
FileList = os.listdir(dataFolder)


for File in FileList:
    if not File.endswith('.xlsx'):
        print ("skipping file named", File)
        continue
    df = pd.read_excel(File, header = None)
    #extract information from the header in the raw data from Noldus
    NumberofHeaderRows = int(df.iloc[0,1])
    Header = df.iloc[31:NumberofHeaderRows - 3,:]
    SubjectFilters = ["Subject", "Mouse","subject", "mouse", "Animal"]
    Subject =  Header[Header[0].isin(SubjectFilters)].iloc[0,1]
#    GenotypeFilters = ["Genotype", "Group", "genotype", "group", "genotype/group"]
#    Genotype =  Header[Header[0].isin(GenotypeFilters)].iloc[0,1]
    DrugFilters = ["Drug", "drug"]
    Drug =  Header[Header[0].isin(DrugFilters)].iloc[0,1]
#    InjFilters = ["Injection time", "injection time", "InjectionTime"]
#    Injection =  Header[Header[0].isin(InjFilters)].iloc[0,1]
    DietFilters = ["Diet", "diet"]
    Diet =  Header[Header[0].isin(DietFilters)].iloc[0,1]
    #GenderFilters = ["Sex", "Gender", "sex", "gender"]
    #Gender =  Header[Header[0].isin(GenderFilters)].iloc[0,1]
    TimestampFilters = ["Timestamp", "timestamp", "<User-defined 1>"]
    Timestamp =  Header[Header[0].isin(TimestampFilters)].iloc[0,1]
    NotesFilters = ["Notes", "notes", "Note"]
    Notes =  Header[Header[0].isin(NotesFilters)].iloc[0,1]
    #print ("Subject:",Subject, "Genotype:", Genotype,"Gender:", Gender,"Timestamp:", Timestamp, "Notes", Notes)
    ColumnNames = df.iloc[[NumberofHeaderRows - 2],:].values.tolist()
    df.columns = ColumnNames
   
    
    #set up a data block to operate on and clean the data
    if df['Trial time'].iloc[-1] > 601:
        DataBlock = df.loc[NumberofHeaderRows +1801:,'Trial time':'Results']
    else:
        DataBlock = df.loc[NumberofHeaderRows +1:,'Trial time':'Results']
    firstRow = DataBlock[:1]
    lastRow = DataBlock[-1:]
    
    firstRow[firstRow == '-'] = 0
    lastRow[lastRow == '-'] = 0
    
    DataBlock[:1] = firstRow
    DataBlock[-1:] = lastRow
     #replace missing data (-) with NaN, then interpolate
    DataBlock.replace('-',np.NaN,inplace = True)
    DataBlock.interpolate(method = 'values', axis = 0, inplace = True)
    isMovingData = DataBlock.loc[:,'Movement(Moving / Center-point)']    
    isMovingData[isMovingData >= 0.5] = 1
    isMovingData[isMovingData < 0.5] = 0
    isInOpen = DataBlock.loc[:,'In zone']
    isInOpen[isInOpen >= 0.5] = 1
    isInOpen[isInOpen < 0.5] = 0
    Velocity = DataBlock.loc [:, 'Velocity']
    
    PercentTimeinOpen = sum (isInOpen[isInOpen == 1])/len(isInOpen)*100
    AvgVelocity = sum (Velocity)/len(Velocity)
    
    
    #slice the data based on the starting and ending rows for movement
    isMovingData.iat[0] = 0
    isMovingData.iat[-1] = 0
    movingTrans = isMovingData.diff()
    startingPoints = movingTrans[movingTrans == 1].index.tolist()
    endingPoints = movingTrans[movingTrans == -1].index.tolist()
    OpenTrue = 0
    inOpenCutOff = 15
    FramesOpenMoving = 0
    MovementBlocks = 0
    TotalMoveDuration = 0
    TotalVelocity = 0
    TotalFramesMoving = 0
   
    
    #loop over the movement data and identify occurances of mouse entering open
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
        isOpenSpan = isInOpen[currentStart:currentEnd]
        numOpen = len(isOpenSpan[isOpenSpan == 1])
        if numOpen > inOpenCutOff:
            OpenTrue +=1
            FramesOpenMoving = FramesOpenMoving + numOpen
    PercentTimeinOpenMoving = FramesOpenMoving/len(isInOpen)*100
    AvgMoveDuration = TotalMoveDuration/MovementBlocks
    AvgVelocityMoving = TotalVelocity/TotalFramesMoving

    
    
    Dataz = [Subject, Drug, Diet, Timestamp, Notes, PercentTimeinOpen, PercentTimeinOpenMoving, MovementBlocks, AvgMoveDuration, TotalMoveDuration, AvgVelocityMoving, AvgVelocity, OpenTrue]
    myDataList.append(Dataz)
    print (myDataList)
    
resultsDF = pd.DataFrame(data = myDataList, columns=resultsColumns)
os.chdir("/Users/leblanckh/data/DREADD_ZeroMaze_RawData/OutputFiles")
writer = pd.ExcelWriter('DREADD_Zero_Maze_Movement_Analysis.xlsx')
resultsDF.to_excel(writer,'Sheet1')
writer.save()


