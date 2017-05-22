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
resultsColumns = ["Subject", "Group", "LED Power", "Timestamp", "Notes", "Time in open, LED off (%)",\
"Time in open, LED on (%)",  "Time in open while moving and LED off(%)", "Time in open while moving and LED on(%)",\
"Number of movements, LED off", "Number of movements, LED on","Average duration of movements, LED off (s)",\
 "Average duration of movements, LED on (s)", "Total time moving, LED off (s)", "Total time moving, LED on (s)",\
 "Speed while moving, LED off (cm/s)", "Speed while moving, LED on (cm/s)","Average velocity, LED off (cm/s)",\
"Average velocity, LED on (cm/s)", "Number of movements in open, LED off", "Number of movements in open, LED on"]
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
    DataBlock = df.loc[NumberofHeaderRows +1801:,'Trial time':'LED OFF']
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
    LEDon = DataBlock.loc[:,'LED ON']
    LEDon[LEDon >= 0.5] = 1
    LEDon[LEDon < 0.5] = 0
    LEDoff = DataBlock.loc[:,'LED OFF']
    Velocity = DataBlock.loc [:, 'Velocity']
    
    AvgVelocity = sum (Velocity)/len(Velocity)
    
    LEDonIdx = LEDon[LEDon ==1].index.tolist()
    LEDon_Filtered_InOpen = isInOpen[LEDonIdx]
    LEDon_Filtered_Velocity = Velocity[LEDonIdx]
    PercentTimeinOpenLEDon = sum (LEDon_Filtered_InOpen[0:])/len(LEDonIdx)*100
    AvgVelocityLEDon = sum(LEDon_Filtered_Velocity)/len(LEDonIdx)
    LEDoffIdx = LEDoff[LEDoff ==1].index.tolist()
    LEDoff_Filtered_InOpen = isInOpen[LEDoffIdx]
    LEDoff_Filtered_Velocity = Velocity[LEDoffIdx]
    PercentTimeinOpenLEDoff = sum (LEDoff_Filtered_InOpen[0:])/len(LEDoffIdx)*100
    AvgVelocityLEDoff = sum(LEDoff_Filtered_Velocity)/len(LEDoffIdx)
    
    
    #slice the data based on the starting and ending rows for movement
    isMovingData.iat[0] = 0
    isMovingData.iat[-1] = 0
    
    movingTrans = isMovingData.diff()
    
    MovingstartingPoints = movingTrans[movingTrans == 1].index.tolist()
    MovingendingPoints = movingTrans[movingTrans == -1].index.tolist()
    
    OpenLEDONTrue = 0
    OpenLEDOFFTrue = 0
    inOpenCutOff = 15
    LEDCutOff = 15
    FramesOpenLEDONMoving = 0
    FramesOpenLEDOFFMoving = 0
    MovementBlocksON = 0
    TotalMoveDurationON = 0
    TotalVelocityON = 0
    TotalFramesMovingON = 0
    MovementBlocksOFF = 0
    TotalMoveDurationOFF = 0
    TotalVelocityOFF = 0
    TotalFramesMovingOFF = 0

    
    #loop over the movement data and identify occurances of mouse entering open
    if len(MovingstartingPoints) != len(MovingendingPoints):
        print ("Uneven start and end pairs. There are", (len(MovingstartingPoints)), "starting points and" , (len(MovingendingPoints)) ,"endingPoints" )
    for i in range(len(MovingstartingPoints)):
        currentStart = MovingstartingPoints[i]+1
        currentEnd = MovingendingPoints[i]
        LEDonSpan = LEDon.loc[currentStart:currentEnd]
        numON = len(LEDonSpan[LEDonSpan ==1])
        numOFF = len(LEDonSpan[LEDonSpan ==0])
        if numON > LEDCutOff:
            LEDonTrans = LEDonSpan.diff()
            LEDonStartingPoints = LEDonTrans[LEDonTrans == 1].index.tolist()
            LEDonEndingPoints = LEDonTrans[LEDonTrans == -1].index.tolist()
            LEDonStartShifted = [x+1 for x in LEDonStartingPoints]
            LEDonEndShifted = [x+1 for x in LEDonEndingPoints]
            if LEDonSpan.iloc[0]==1:
                LEDonStartShifted.insert(0,currentStart)
            if LEDonSpan.iloc[-1]==1:
                LEDonEndShifted.append(currentEnd)
            for i in range(len(LEDonStartShifted)):
                currentLEDonStart = LEDonStartShifted[i]
                currentLEDonEnd = LEDonEndShifted[i]
                MovementBlocksON +=1
                MovementDurationON = (currentLEDonEnd-currentLEDonStart)/30
                TotalMoveDurationON = TotalMoveDurationON + MovementDurationON
                TotalFramesMovingON = TotalFramesMovingON + numON
                MoveVelocityON = sum (Velocity.loc[currentLEDonStart:currentLEDonEnd])
                TotalVelocityON = TotalVelocityON + MoveVelocityON
                isOpenSpan = isInOpen.loc[currentLEDonStart:currentLEDonEnd]
                numOpen = len(isOpenSpan[isOpenSpan == 1])
                if numOpen > inOpenCutOff:
                    OpenLEDONTrue +=1
                    FramesOpenLEDONMoving = FramesOpenLEDONMoving + numOpen
        if numOFF > LEDCutOff:
            LEDoffTrans = LEDonSpan.diff()
            LEDoffStartingPoints = LEDoffTrans[LEDoffTrans == -1].index.tolist()
            LEDoffEndingPoints = LEDoffTrans[LEDoffTrans == 1].index.tolist()
            LEDoffStartShifted = [x+1 for x in LEDoffStartingPoints]
            LEDoffEndShifted = [x+1 for x in LEDoffEndingPoints]
            if LEDonSpan.iloc[0]==0:
                LEDoffStartShifted.insert(0,currentStart)
            if LEDonSpan.iloc[-1]==0:
                LEDoffEndShifted.append(currentEnd)
            for i in range(len(LEDoffStartShifted)):
                currentLEDoffStart = LEDoffStartShifted[i]
                currentLEDoffEnd = LEDoffEndShifted[i]
#                print("moveStart", currentStart, "moveEnd", currentEnd, "LEDoffstart", currentLEDoffStart, 'LEDoffend', currentLEDoffEnd)
                MovementBlocksOFF +=1
                MovementDurationOFF = (currentLEDoffEnd-currentLEDoffStart)/30
                TotalMoveDurationOFF = TotalMoveDurationOFF + MovementDurationOFF
                TotalFramesMovingOFF = TotalFramesMovingOFF + numOFF
                MoveVelocityOFF = sum (Velocity.loc[currentLEDoffStart:currentLEDoffEnd])
                TotalVelocityOFF = TotalVelocityOFF + MoveVelocityOFF
                isOpenSpan = isInOpen.loc[currentLEDoffStart:currentLEDoffEnd]
                numOpen = len(isOpenSpan[isOpenSpan == 1])
                if numOpen > inOpenCutOff:
                    OpenLEDOFFTrue +=1
                    FramesOpenLEDOFFMoving = FramesOpenLEDOFFMoving + numOpen
#                    print ("isOpenSpan", isOpenSpan, "currentStart", currentLEDoffStart, "currentEnd", currentLEDoffEnd, "Open entries", OpenLEDOFFTrue)

       
    PercentTimeinOpenMovingON = FramesOpenLEDONMoving/len(isInOpen)*100
    PercentTimeinOpenMovingOFF = FramesOpenLEDOFFMoving/len(isInOpen)*100
    AvgMoveDurationON = TotalMoveDurationON/MovementBlocksON
    AvgMoveDurationOFF = TotalMoveDurationOFF/MovementBlocksOFF
    AvgVelocityMovingON = TotalVelocityON/TotalFramesMovingON
    AvgVelocityMovingOFF = TotalVelocityOFF/TotalFramesMovingOFF
#    print("AvgMoveDurationON, OFF", AvgMoveDurationON, AvgMoveDurationOFF, \
#    "AvgVelocityON,OFF",AvgVelocityMovingON, AvgVelocityMovingOFF, \
#    "MoveON,OFF", MovementBlocksON, MovementBlocksOFF, \
#    "PercentTimeinOpenMovingON,OFF", PercentTimeinOpenMovingON, PercentTimeinOpenMovingOFF,\
#    "NumberofMoveinOpenON,OFF", OpenLEDONTrue, OpenLEDOFFTrue)
#    
#    
    Dataz = [Subject,Genotype,LEDPower, Timestamp, Notes, PercentTimeinOpenLEDoff,\
    PercentTimeinOpenLEDon, PercentTimeinOpenMovingOFF, PercentTimeinOpenMovingON,\
    MovementBlocksOFF, MovementBlocksON, AvgMoveDurationOFF, AvgMoveDurationON,\
    TotalMoveDurationOFF, TotalMoveDurationON, AvgVelocityMovingOFF, AvgVelocityMovingON,\
    AvgVelocityLEDoff, AvgVelocityLEDon, OpenLEDOFFTrue, OpenLEDONTrue]
    myDataList.append(Dataz)
    print (myDataList)
    
resultsDF = pd.DataFrame(data = myDataList, columns=resultsColumns)
os.chdir("/Users/leblanckh/data/Opto_ZeroMaze_RawData/OutputFiles")
writer = pd.ExcelWriter('Opto_Zero_Maze_Movement_Analysis.xlsx')
resultsDF.to_excel(writer,'Sheet1')
writer.save()


