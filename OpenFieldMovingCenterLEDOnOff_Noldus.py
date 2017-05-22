# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 21:28:53 2017

@author: james.leblanc
"""

import pandas as pd
import numpy as np
import os
import sys
dataFolder = "/Users/leblanckh/data/Opto_OpenField_RawData"
resultsColumns = ["Subject", "Group", "LED Power", "Timestamp", "Notes", "Time in Center, LED off (%)",\
"Time in Center, LED on (%)",  "Time in Center while moving and LED off(%)", "Time in Center while moving and LED on(%)",\
"Number of movements, LED off", "Number of movements, LED on","Average duration of movements, LED off (s)",\
 "Average duration of movements, LED on (s)", "Total time moving, LED off (s)", "Total time moving, LED on (s)",\
 "Speed while moving, LED off (cm/s)", "Speed while moving, LED on (cm/s)","Average velocity, LED off (cm/s)",\
"Average velocity, LED on (cm/s)", "Number of movements in Center, LED off", "Number of movements in Center, LED on"]
myDataList = []
os.chdir(dataFolder)
FileList = os.listdir(dataFolder)

ONE_MINUTE_AS_FRAMES = 1798
TWO_MINUTES_AS_FRAMES = 2 * ONE_MINUTE_AS_FRAMES

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
        DataBlock = df.loc[NumberofHeaderRows +8992:,'Trial time':'LED OFF']
    
    LEDon = DataBlock.loc[:,'LED ON']
    
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
    
    #slice the data based on the starting and ending rows for movement
    isMovingData.iat[0] = 0
    isMovingData.iat[-1] = 0
    movingTrans = isMovingData.diff()
    MovingstartingPoints = movingTrans[movingTrans == 1].index.tolist()
    MovingendingPoints = movingTrans[movingTrans == -1].index.tolist()
    print("before fix: Move start,end", MovingstartingPoints, MovingendingPoints)
    
# for trials with 2 minute LED off periods, discard the 2nd minute of data and reassign DataBlock variables
    if df['Trial time'].iloc[-1] >601:
        TotalNumberDataRows = len(LEDon)
        FinalIndexLabel = DataBlock[-1:].index.tolist()[0]
        LEDon.iat[0] = 0
        LEDon.iat[-1] = 0
        LEDtrans = LEDon.diff()
        LEDoffStart = LEDtrans[LEDtrans == -1].index.tolist()
        LEDoffMin2Start = [x+ONE_MINUTE_AS_FRAMES for x in LEDoffStart]
        LEDoffMin2End = [x+TWO_MINUTES_AS_FRAMES for x in LEDoffStart]
        AllMinutes2drop = set()

        for i in range(len(LEDoffMin2Start)):
            curStart = LEDoffMin2Start[i]+1
            curEnd = LEDoffMin2End[i]+1
            print("current start,end", curStart,curEnd)

            #protect against End of File (EoF)
            if curStart >= FinalIndexLabel:
                break;
            else:
                if curEnd >= FinalIndexLabel:
                    curEnd = FinalIndexLabel - 1 #if EOF need to adjust down 1
                    
            for i in range(len(MovingstartingPoints)):
                MoveStart = MovingstartingPoints[i]+1
                MoveEnd = MovingendingPoints[i]
                if MoveEnd in range (curStart,curEnd):
                    print("For MoveEnd", MoveEnd, "before: move column at curStart", DataBlock.loc[curStart-1,'Movement(Moving / Center-point)'], curStart-1)
                    DataBlock.loc[curStart-1,'Movement(Moving / Center-point)'] = 0
                    print("For MoveEnd", MoveEnd, "after:move column at curStart", DataBlock.loc[curStart-1,'Movement(Moving / Center-point)'], curStart-1)
                if MoveStart in range(curStart,curEnd):
                    print("For MoveStart", MoveStart, "before: move column at curEnd", DataBlock.loc[curEnd+1,'Movement(Moving / Center-point)'], curEnd+1)
                    DataBlock.loc[curEnd+1,'Movement(Moving / Center-point)'] = 0
                    print("For MoveStart", MoveStart, "after: move column at curEnd", DataBlock.loc[curEnd+1,'Movement(Moving / Center-point)'], curEnd+1)

                

            print("JJ FinalIdxL is ", FinalIndexLabel)
            #create a set that includes all index labels for the current minute
            if curEnd<=FinalIndexLabel:
                Minute2dropSpan = set(range(curStart, curEnd+1))
                print(" XXXXXXX ")
                print("Minute2dropSpan max val is ", max(Minute2dropSpan))
            else:
                Minute2dropSpan = set(range(curStart, curEnd))
                print(" XXXXXXX ")
                print("Minute2dropSpan max val is ", max(Minute2dropSpan))
            #perform set Union to accumulate the superset that contains all minutes to be dropped
            AllMinutes2drop = AllMinutes2drop | Minute2dropSpan  

        # Outside the loop create a set of index labels that should be kept
        idx_2_keep = set(DataBlock.index.tolist()) - AllMinutes2drop
        idx_2_keepAsList = list(idx_2_keep)
        idx_2_keepAsList.sort()
#        print(" *****  Before all of that nonsense  ***** ")
##        print("idx_2_keepAsList is ", idx_2_keepAsList)
#        print("idx_2_keepAsList[0] is ", idx_2_keepAsList[0], " and idx_2_keepAsList[-1] is ", idx_2_keepAsList[-1])
#        print(" *****  Before the suplex  ***** ")        
#        print("DataBlock.iat[0].index is ", DataBlock.index.tolist()[0])
#        print("DataBlock.iat[-1].index is ", DataBlock.index.tolist()[-1])
#        print("Length of Datablock", len(LEDon))        
        DataBlock = DataBlock.loc[idx_2_keepAsList]
        print(" *****  After the suplex  ***** ")
        print("DB moving around start of interest (12625) is", DataBlock.loc[12625:14440,'Movement(Moving / Center-point)'])
#        print("DataBlock.iat[0].index is ", DataBlock.index.tolist()[0])
#        print("DataBlock.iat[-1].index is ", DataBlock.index.tolist()[-1])
#        print("length of DataBlock", len(LEDon))
        isMovingData = DataBlock.loc[:,'Movement(Moving / Center-point)']
        isMovingData.iat[-1] = 0
        movingTrans = isMovingData.diff()
        MovingstartingPoints = movingTrans[movingTrans == 1].index.tolist()
        MovingendingPoints = movingTrans[movingTrans == -1].index.tolist()
        print("after fix: Move start,end", MovingstartingPoints, MovingendingPoints)


    isMovingData = DataBlock.loc[:,'Movement(Moving / Center-point)']
    isInCenter = DataBlock.loc[:,'In zone']
    LEDon = DataBlock.loc[:,'LED ON']
    TotalNumberDataRows = len(LEDon)
    LEDoff = DataBlock.loc[:,'LED OFF']
    Velocity = DataBlock.loc [:, 'Velocity']
    PercentTimeinCenter = sum (isInCenter[isInCenter == 1])/len(isInCenter)*100
    AvgVelocity = sum (Velocity)/len(Velocity)
    
    LEDonIdx = LEDon[LEDon ==1].index.tolist()
    LEDon_Filtered_InCenter = isInCenter[LEDonIdx]
    LEDon_Filtered_Velocity = Velocity[LEDonIdx]
    PercentTimeinCenterLEDon = sum (LEDon_Filtered_InCenter[0:])/len(LEDonIdx)*100
    AvgVelocityLEDon = sum(LEDon_Filtered_Velocity)/len(LEDonIdx)
    LEDoffIdx = LEDoff[LEDoff ==1].index.tolist()
    LEDoff_Filtered_InCenter = isInCenter[LEDoffIdx]
    LEDoff_Filtered_Velocity = Velocity[LEDoffIdx]
    PercentTimeinCenterLEDoff = sum (LEDoff_Filtered_InCenter[0:])/len(LEDoffIdx)*100
    AvgVelocityLEDoff = sum(LEDoff_Filtered_Velocity)/len(LEDoffIdx)
    
    
    CenterLEDONTrue = 0
    CenterLEDOFFTrue = 0
    inCenterCutOff = 15
    LEDCutOff = 15
    FramesCenterLEDONMoving = 0
    FramesCenterLEDOFFMoving = 0
    MovementBlocksON = 0
    TotalMoveDurationON = 0
    TotalVelocityON = 0
    TotalFramesMovingON = 0
    MovementBlocksOFF = 0
    TotalMoveDurationOFF = 0
    TotalVelocityOFF = 0
    TotalFramesMovingOFF = 0

    
    #loop over the movement data and identify occurances of mouse entering Center
    if len(MovingstartingPoints) != len(MovingendingPoints):
        print ("Uneven start and end pairs. There are", (len(MovingstartingPoints)), "starting points and" , (len(MovingendingPoints)) ,"endingPoints" )
    for i in range(len(MovingstartingPoints)):
        currentStart = MovingstartingPoints[i]+1
        currentEnd = MovingendingPoints[i]
        print("current start, end", currentStart, currentEnd)
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
                isCenterSpan = isInCenter.loc[currentLEDonStart:currentLEDonEnd]
                numCenter = len(isCenterSpan[isCenterSpan == 1])
                if numCenter > inCenterCutOff:
                    CenterLEDONTrue +=1
                    FramesCenterLEDONMoving = FramesCenterLEDONMoving + numCenter
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
                isCenterSpan = isInCenter.loc[currentLEDoffStart:currentLEDoffEnd]
                numCenter = len(isCenterSpan[isCenterSpan == 1])
                if numCenter > inCenterCutOff:
                    CenterLEDOFFTrue +=1
                    FramesCenterLEDOFFMoving = FramesCenterLEDOFFMoving + numCenter
#                    print ("isCenterSpan", isCenterSpan, "currentStart", currentLEDoffStart, "currentEnd", currentLEDoffEnd, "Center entries", CenterLEDOFFTrue)

       
    PercentTimeinCenterMovingON = FramesCenterLEDONMoving/len(isInCenter)*100
    PercentTimeinCenterMovingOFF = FramesCenterLEDOFFMoving/len(isInCenter)*100
    AvgMoveDurationON = TotalMoveDurationON/MovementBlocksON
    AvgMoveDurationOFF = TotalMoveDurationOFF/MovementBlocksOFF
    AvgVelocityMovingON = TotalVelocityON/TotalFramesMovingON
    AvgVelocityMovingOFF = TotalVelocityOFF/TotalFramesMovingOFF
    
    
    Dataz = [Subject,Genotype,LEDPower, Timestamp, Notes, PercentTimeinCenterLEDoff,\
    PercentTimeinCenterLEDon, PercentTimeinCenterMovingOFF, PercentTimeinCenterMovingON,\
    MovementBlocksOFF, MovementBlocksON, AvgMoveDurationOFF, AvgMoveDurationON,\
    TotalMoveDurationOFF, TotalMoveDurationON, AvgVelocityMovingOFF, AvgVelocityMovingON,\
    AvgVelocityLEDoff, AvgVelocityLEDon, CenterLEDOFFTrue, CenterLEDONTrue]
    myDataList.append(Dataz)
    print (myDataList)
    
resultsDF = pd.DataFrame(data = myDataList, columns=resultsColumns)
os.chdir("/Users/leblanckh/data/Opto_OpenField_RawData/OutputFiles")
writer = pd.ExcelWriter('Opto_OpenField_Movement_Analysis.xlsx')
resultsDF.to_excel(writer,'Sheet1')
writer.save()





