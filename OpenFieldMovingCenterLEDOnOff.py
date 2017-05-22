# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 21:28:53 2017

@author: leblanckh
"""

import pandas as pd
import numpy as np
import os
import sys
dataFolder = "/Users/leblanckh/data/Opto_ZeroMaze_RawData"
resultsColumns = ["Subject", "Group", "LED Power", "Timestamp", "Notes", "Time in Center (%)", "Time in Center while moving and LED ON(%)", "Time in Center while moving and LED OFF(%)","Number of movements", "Average duration of movements (s)", "Total time moving (s)", "Speed while moving (cm/s)", "Average velocity (cm/s)", "Number of movements in Center, LED ON", "Number of movements in Center, LED OFF"]
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
    GenotypeFilters = ["Genotype", "Group", "genotype", "group", "genotype/group"]
    Genotype =  Header[Header[0].isin(GenotypeFilters)].iloc[0,1]
    LEDFilters = ["LED Stimulation", "LED power", "LED power (mW)"]
    LEDPower = Header[Header[0].isin(LEDFilters)].iloc[0,1]
    TimestampFilters = ["Timestamp", "timestamp", "<User-defined 1>"]
    Timestamp =  Header[Header[0].isin(TimestampFilters)].iloc[0,1]
    NotesFilters = ["Notes", "notes", "Note"]
    Notes =  Header[Header[0].isin(NotesFilters)].iloc[0,1]
    #print ("Subject:",Subject, "Genotype:", Genotype,"LEDPower:", LEDPower, "Timestamp:", Timestamp, "Notes", Notes)
    ColumnNames = df.iloc[[NumberofHeaderRows - 2],:].values.tolist()
    df.columns = ColumnNames
   
    
    #set up a data block to operate on and clean the data
    if df['Trial time'].iloc[-1] <601:
        DataBlock = df.loc[NumberofHeaderRows +1:,'Trial time':'LED OFF']
    else:
        DataBlock = df.loc[NumberofHeaderRows +10801:,'Trial time':'LED OFF']
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
    isInCenter = DataBlock.loc[:,'In zone']
    LEDon = DataBlock.loc[:,'LED ON']
    LEDon[LEDon >= 0.5] = 1
    LEDon[LEDon < 0.5] = 0
    LEDoff = DataBlock.loc[:,'LED OFF']
    Velocity = DataBlock.loc [:, 'Velocity']
    
    PercentTimeinCenter = sum (isInCenter[isInCenter == 1])/len(isInCenter)*100
    AvgVelocity = sum (Velocity)/len(Velocity)
    
    #slice the data based on the starting and ending rows for movement
    isMovingData.iat[0] = 0
    isMovingData.iat[-1] = 0
    movingTrans = isMovingData.diff()
    startingPoints = movingTrans[movingTrans == 1].index.tolist()
    endingPoints = movingTrans[movingTrans == -1].index.tolist()
    CenterLEDONTrue = 0
    CenterLEDOFFTrue = 0
    inCenterCutOff = 15
    FramesCenterLEDONMoving = 0
    FramesCenterLEDOFFMoving = 0
    MovementBlocks = 0
    TotalMoveDuration = 0
    TotalVelocity = 0
    TotalFramesMoving = 0

    
    #loop over the movement data and identify occurances of mouse entering Center
    if len(startingPoints) != len(endingPoints):
        print ("Uneven start and end pairs. There are", (len(startingPoints)), "starting points and" , (len(endingPoints)) ,"endingPoints" )
    for i in range(len(startingPoints)):
        currentStart = startingPoints[i]
        currentEnd = endingPoints[i]-1
        LEDonMove = LEDon.loc[currentStart]
        MovementBlocks += 1
        MovementDuration = (currentEnd-currentStart)/30
        TotalMoveDuration = TotalMoveDuration + MovementDuration
        FramesMoving = (currentEnd-currentStart)
        TotalFramesMoving = TotalFramesMoving + FramesMoving
        MoveVelocity = sum (Velocity.loc[currentStart:currentEnd])
        TotalVelocity = TotalVelocity + MoveVelocity
#        print ("LEDonSpan", LEDonSpan)
#        LEDonEpochs = LEDonSpan[LEDonSpan == 1].index.tolist()
#        print ("LEDonEpochs", LEDonEpochs)
#        isCenterLEDon = isInCenter[LEDonEpochs]
#        print ("isCenterLEDon", isCenterLEDon)
#        LEDoffSpan = LEDoff[currentStart:currentEnd]
        isCenterSpan = isInCenter[currentStart:currentEnd]
        numCenter = len(isCenterSpan[isCenterSpan == 1])
        if numCenter > inCenterCutOff:
            if LEDonMove == 1:
                CenterLEDONTrue +=1
                FramesCenterLEDONMoving = FramesCenterLEDONMoving + numCenter
            else:
                CenterLEDOFFTrue +=1
                FramesCenterLEDOFFMoving = FramesCenterLEDOFFMoving + numCenter
    PercentTimeinCenterMovingON = FramesCenterLEDONMoving/len(isInCenter)*100
    PercentTimeinCenterMovingOFF = FramesCenterLEDOFFMoving/len(isInCenter)*100
    AvgMoveDuration = TotalMoveDuration/MovementBlocks
    AvgVelocityMoving = TotalVelocity/TotalFramesMoving
    
    
    Dataz = [Subject,Genotype,LEDPower, Timestamp, Notes, PercentTimeinCenter, PercentTimeinCenterMovingON, PercentTimeinCenterMovingOFF, MovementBlocks, AvgMoveDuration, TotalMoveDuration, AvgVelocityMoving, AvgVelocity, CenterLEDONTrue, CenterLEDOFFTrue]
    myDataList.append(Dataz)
    print (myDataList)
    
resultsDF = pd.DataFrame(data = myDataList, columns=resultsColumns)
os.chdir("/Users/leblanckh/data/Opto_OpenField_RawData/OutputFiles")
writer = pd.ExcelWriter('Opto_OpenField_Movement_Analysis.xlsx')
resultsDF.to_excel(writer,'Sheet1')
writer.save()


