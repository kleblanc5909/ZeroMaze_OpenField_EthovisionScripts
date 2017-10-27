# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 21:28:53 2017

@author: leblanckh
"""

import pandas as pd
import numpy as np
import os
import sys
dataFolder = "/Users/leblanckh/data/Opto_OpenField_RawData"
resultsColumns = ["Subject", "Group", "LED Power", "Timestamp", "Notes", "Time in Center, LED on (%)",\
"Time in Center, LED off (%)",  "Time in Center while moving and LED on(%)", "Time in Center while moving and LED off(%)",\
"Time in Surround while moving and LED on(%)", "Time in Surround while moving and LED off(%)", \
"Time moving in Center only, LED on (%)", "Time moving in Center only, LED off (%)",\
"Time moving in Surround only, LED on (%)", "Time moving in Surround only, LED off (%)",\
"Number of movements, LED on", "Number of movements, LED off","Average duration of movements, LED on (s)",\
 "Average duration of movements, LED off (s)", "Total time moving, LED on (s)", "Total time moving, LED off (s)",\
 "Speed while moving, LED on (cm/s)", "Speed while moving, LED off (cm/s)",\
 "Average Velocity for Center Only Movements, LED on (cm/s)","Average Velocity for Center Only Movements, LED off (cm/s)", \
 "Average Velocity for Surround Only Movements, LED on (cm/s)", "Average Velocity for Surround Only Movements, LED off (cm/s)", \
 "Average Velocity for Center Movements, LED on (cm/s)","Average Velocity for Center Movements, LED off (cm/s)", \
 "Average Velocity for Surround Movements, LED on (cm/s)","Average Velocity for Surround Movements, LED off (cm/s)", \
 "Average velocity, LED on (cm/s)","Average velocity, LED off (cm/s)", \
 "Number of movements in Center, LED on", "Number of movements in Center, LED off", \
 "Number of movements in Surround, LED on","Number of movements in Surround, LED off",\
 "Number of Center-only movements, LED on", "Number of Center-only movements, LED off", \
 "Number of Surround-only movements, LED on", "Number of Surround-only movements, LED off"]
 
myDataList = []
os.chdir(dataFolder)
FileList = os.listdir(dataFolder)

ONE_MINUTE_AS_FRAMES = 1798
TWO_MINUTES_AS_FRAMES = 2 * ONE_MINUTE_AS_FRAMES
NO_BASELINE_10MIN_SESSION_TRIAL_LENGTH = 601

def Extract_Header_Info(dataframe):
    """
    Identifies the length of the header,Extracts independent variable 
    information from the header, and assigns the column names to the row that
    contains the column information in the header
    
    inputs: dataframe
    Returns: NumberofHeaderRows, Subject, Drug, Diet, Timestamp, Notes, 
    """
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
    
    return NumberofHeaderRows,Subject,Genotype,LEDPower,Timestamp,Notes

def create_DataBlock(dataframe, FirstCutOffTime, shortTimeIdx, longTimeIdx):
    """
    Creates a dataFrame, theData, that contains the actual raw data values.
    It also replaces missing values in the movement column of data with NaN, then interpolates
    returns: theData, movement column data (isMovingData), LED on column (LEDon)
    """
    finalColumnName = dataframe.columns[-1]
    initialColumnName = dataframe.columns[0]
    
    if dataframe['Trial time'].iloc[-1] <FirstCutOffTime:
        theData = dataframe.loc[shortTimeIdx:,initialColumnName:finalColumnName] #NumHdr + 1
    else:
        theData = dataframe.loc[longTimeIdx:,initialColumnName:finalColumnName] #NumHdr + 8992
    
    theData[:1].replace('-',0,inplace = True)
    theData[-1:].replace('-',0,inplace = True)
     #replace missing data (-) with NaN, then interpolate
    theData.replace('-',np.NaN,inplace = True)
    theData.interpolate(method = 'values', axis = 0, inplace = True)
    theData.loc[:,'Movement(Moving / Center-point)'] = Binary_Data_Interpolation_CleanUp(theData.loc[:,'Movement(Moving / Center-point)']) 
    theData.loc[:,'In zone'] = Binary_Data_Interpolation_CleanUp(theData.loc[:,'In zone']) 
    theData.loc[:,'LED ON'] = Binary_Data_Interpolation_CleanUp(theData.loc[:,'LED ON']) 
    
    return theData
    
def Binary_Data_Interpolation_CleanUp(binaryDataColumn):
    """
    After interpolation, binary data columns need decimal point values converted to 0s or 1s.
    This function sets all values 0.5 and greater to 1, and all values less than 0.5 to 0.
    Inputs: the binary data column with decimal values
    Returns: the binary data column with 0s and 1s only
    """
    
    binaryDataColumn[binaryDataColumn >= 0.5] = 1
    binaryDataColumn[binaryDataColumn < 0.5] = 0
    
    return binaryDataColumn

def Binary_Data_Transition_Point_Finder(binaryDataColumn):
    """
    Uses a diff method to find when 0 changes to 1 and vice versa in binary data column
    Defines these transitions as either the starting point or ending point of an action or state, respectively
    and sends that index value to a list
    Also sets the first and last value of the column to 0 for edge case handling
    
    input: binaryDataColumn
    Returns:  a list of binaryDataStartPoints and binaryDataEndPoints
    """
    binaryDataColumn.iat[0] = 0
    binaryDataColumn.iat[-1] = 0
    binaryDataTrans = binaryDataColumn.diff()
    binaryDataStartingPoints = binaryDataTrans[binaryDataTrans == 1].index.tolist()
    binaryDataEndingPoints = binaryDataTrans[binaryDataTrans == -1].index.tolist()
    
    return binaryDataStartingPoints, binaryDataEndingPoints
    
def Remove_Minute2_LED_OFF (theData, MoveStart, MoveEnd):
    """
    for trials with 2 minute LED off periods, discard the 2nd minute of data and reassign DataBlock variables
    
    inputs: dataframe and the starting and ending points of the moving column
    returns: adjusted start and end of movements
    """
    #TotalNumberDataRows = len(LEDon)
    FinalIndexLabel = theData[-1:].index.tolist()[0]
    LEDonStart,LEDoffStart = Binary_Data_Transition_Point_Finder(LEDon)
    LEDoffMin2Start = [x+ONE_MINUTE_AS_FRAMES for x in LEDoffStart]
    LEDoffMin2End = [x+TWO_MINUTES_AS_FRAMES for x in LEDoffStart]
    AllMinutes2drop = set()

    for i in range(len(LEDoffMin2Start)):
        curStart = LEDoffMin2Start[i]+1
        curEnd = LEDoffMin2End[i]+1


        #protect against End of File (EoF)
        if curStart >= FinalIndexLabel:
            break;
        else:
            if curEnd >= FinalIndexLabel:
                curEnd = FinalIndexLabel - 1 #if EOF need to adjust down 1
                
        for i in range(len(MoveStart)):
            MovingStart = MoveStart[i]+1
            MovingEnd = MoveEnd[i]
            if MovingEnd in range (curStart,curEnd):
                    print("For MoveEnd", MovingEnd, "before: move column at curStart", theData.loc[curStart-1,'Movement(Moving / Center-point)'], curStart-1)
                    theData.loc[curStart-1,'Movement(Moving / Center-point)'] = 0
                    print("For MoveEnd", MovingEnd, "after:move column at curStart", theData.loc[curStart-1,'Movement(Moving / Center-point)'], curStart-1)
            if MovingStart in range(curStart,curEnd):
                    print("For MoveStart", MovingStart, "before: move column at curEnd", theData.loc[curEnd+1,'Movement(Moving / Center-point)'], curEnd+1)
                    theData.loc[curEnd+1,'Movement(Moving / Center-point)'] = 0
                    print("For MoveStart", MovingStart, "after: move column at curEnd", theData.loc[curEnd+1,'Movement(Moving / Center-point)'], curEnd+1)

            

#            print("JJ FinalIdxL is ", FinalIndexLabel)
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
        print("All Minutes to drop", AllMinutes2drop)

    # Outside the loop create a set of index labels that should be kept
    idx_2_keep = set(theData.index.tolist()) - AllMinutes2drop
    idx_2_keepAsList = list(idx_2_keep)
    idx_2_keepAsList.sort()
    print(" *****  Before all of that nonsense  ***** ")
    print("idx_2_keepAsList is ", idx_2_keepAsList)
    print("idx_2_keepAsList[0] is ", idx_2_keepAsList[0], " and idx_2_keepAsList[-1] is ", idx_2_keepAsList[-1])
    print(" *****  Before the suplex  ***** ")        
    print("theData.iat[0].index is ", theData.index.tolist()[0])
    print("theData.iat[-1].index is ", theData.index.tolist()[-1])
    print("Length of Datablock", len(LEDon))        
    theData = theData.loc[idx_2_keepAsList]
    print(" *****  After the suplex  ***** ")
    print("DB moving around start of interest (12625) is", theData.loc[12625:14440,'Movement(Moving / Center-point)'])
    print("DataBlock.iat[0].index is ", theData.index.tolist()[0])
    print("DataBlock.iat[-1].index is ", theData.index.tolist()[-1])
    print("length of DataBlock", len(LEDon))
    isMovingData = theData.loc[:,'Movement(Moving / Center-point)']
    isMovingData.iat[-1] = 0
    movingTrans = isMovingData.diff()
    NewMoveStart = movingTrans[movingTrans == 1].index.tolist()
    NewMoveEnd = movingTrans[movingTrans == -1].index.tolist()
    print("after fix: Move start,end", NewMoveStart, NewMoveEnd)
    
    return NewMoveStart, NewMoveEnd
    
def Movement_Analysis(theData, InZoneColumn,VelocityColumn,NewMoveStart,NewMoveEnd):
    """
    Sets up the initial conditions, loops over the movement data and counts movement blocks,
    duration and velocity,identifies movements into the center, and runs all of the calculations
    
    
    inputs: inZoneColumn, Velocity Column, and the starting and ending points of the moving column
    returns: PercentTimeinCenter, PercentTimeinCenterMoving,MovementBlocks, AvgMoveDuration, 
    TotalMoveDuration, AvgVelocityMoving, AvgVelocity,CenterTrue
    """
    
    #isMovingData = theData.loc[:,'Movement(Moving / Center-point)']
    InZoneColumn = theData.loc[:,'In zone']
    LEDon = theData.loc[:,'LED ON']
    #TotalNumberDataRows = len(LEDon)
    LEDoff = theData.loc[:,'LED OFF']
    VelocityColumn = theData.loc [:, 'Velocity']
    #PercentTimeinCenter = sum (InZoneColumn[InZoneColumn == 1])/len(InZoneColumn)*100
    #AvgVelocity = sum (VelocityColumn)/len(VelocityColumn)
    
    LEDonIdx = LEDon[LEDon ==1].index.tolist()
    LEDon_Filtered_InCenter = InZoneColumn[LEDonIdx]
    LEDon_Filtered_Velocity = VelocityColumn[LEDonIdx]
    PercentTimeinCenterLEDon = sum (LEDon_Filtered_InCenter[0:])/len(LEDonIdx)*100
    AvgVelocityLEDon = sum(LEDon_Filtered_Velocity)/len(LEDonIdx)
    LEDoffIdx = LEDoff[LEDoff ==1].index.tolist()
    LEDoff_Filtered_InCenter = InZoneColumn[LEDoffIdx]
    LEDoff_Filtered_Velocity = VelocityColumn[LEDoffIdx]
    PercentTimeinCenterLEDoff = sum (LEDoff_Filtered_InCenter[0:])/len(LEDoffIdx)*100
    AvgVelocityLEDoff = sum(LEDoff_Filtered_Velocity)/len(LEDoffIdx)
    
    
    CenterLEDONTrue = 0
    CenterLEDOFFTrue = 0
    CenterLEDONOnly = 0
    SurroundLEDONOnly = 0
    SurroundLEDONTrue = 0
    CenterLEDOFFOnly = 0
    SurroundLEDOFFOnly = 0
    SurroundLEDOFFTrue = 0
    inCenterCutOff = 15
    LEDCutOff = 15
    FramesCenterLEDONMoving = 0
    FramesCenterLEDOFFMoving = 0
    FramesCenterLEDONOnly = 0
    FramesSurroundLEDONOnly = 0
    FramesSurroundLEDONMoving = 0
    FramesCenterLEDOFFOnly = 0
    FramesSurroundLEDOFFOnly = 0
    FramesSurroundLEDOFFMoving = 0
    MovementBlocksON = 0
    TotalMoveDurationON = 0
    TotalVelocityON = 0
    TotalVelCenterONOnly = 0
    TotalVelSurroundONOnly = 0
    TotalVelCenterONMove = 0
    TotalVelSurroundONMove = 0
    TotalVelCenterOFFOnly = 0
    TotalVelSurroundOFFOnly = 0
    TotalVelCenterOFFMove = 0
    TotalVelSurroundOFFMove = 0
    TotalFramesMovingON = 0
    MovementBlocksOFF = 0
    TotalMoveDurationOFF = 0
    TotalVelocityOFF = 0
    TotalFramesMovingOFF = 0

    
    #loop over the movement data and identify occurances of mouse entering Center
    if len(NewMoveStart) != len(NewMoveEnd):
        print ("Uneven start and end pairs. There are", (len(NewMoveStart)), "starting points and" , (len(NewMoveEnd)) ,"endingPoints" )
    for i in range(len(NewMoveStart)):
        currentStart = NewMoveStart[i]+1
        currentEnd = NewMoveEnd[i]
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
                print ("current ON start:", currentLEDonStart, "current ON end:", currentLEDonEnd)
                MovementBlocksON +=1
                MovementDurationON = (currentLEDonEnd-currentLEDonStart)/30
                TotalMoveDurationON = TotalMoveDurationON + MovementDurationON
                TotalFramesMovingON = TotalFramesMovingON + numON
                MoveVelocityON = sum (VelocityColumn.loc[currentLEDonStart:currentLEDonEnd])
                TotalVelocityON = TotalVelocityON + MoveVelocityON
                isCenterSpanON = InZoneColumn.loc[currentLEDonStart:currentLEDonEnd]
                numCenterON = len(isCenterSpanON[isCenterSpanON == 1])
                numSurroundON = len(isCenterSpanON[isCenterSpanON == 0])
                inCenterIndexON = isCenterSpanON[isCenterSpanON == 1].index.tolist()
                inSurroundIndexON = isCenterSpanON[isCenterSpanON == 0].index.tolist()
                
                CenterThreshold = 0.98
                SurroundThreshold = 0.02
                PercentCenterON = numCenterON/numON
        
                if PercentCenterON >=CenterThreshold:
                    CenterLEDONOnly +=1
                    FramesCenterLEDONOnly = FramesCenterLEDONOnly + numCenterON
                    VelocityCenterOnly = sum(VelocityColumn.loc[inCenterIndexON])
                    TotalVelCenterONOnly = TotalVelCenterONOnly + VelocityCenterOnly
                    print("Frames Center Only", FramesCenterLEDONOnly)
                if PercentCenterON <=SurroundThreshold:
                    SurroundLEDONOnly +=1
                    FramesSurroundLEDONOnly = FramesSurroundLEDONOnly + numSurroundON
                    VelocitySurroundOnly = sum(VelocityColumn.loc[inSurroundIndexON])
                    TotalVelSurroundONOnly = TotalVelSurroundONOnly + VelocitySurroundOnly
                    print("Frames Surround Only", FramesSurroundLEDONOnly)
                if numCenterON > inCenterCutOff:
                    CenterLEDONTrue +=1
                    FramesCenterLEDONMoving = FramesCenterLEDONMoving + numCenterON
                    VelocityCenterMove = sum(VelocityColumn.loc[inCenterIndexON])
                    TotalVelCenterONMove = TotalVelCenterONMove + VelocityCenterMove
                    print("Frames Center Move", FramesCenterLEDONMoving)
                if numSurroundON > inCenterCutOff:
                    SurroundLEDONTrue +=1
                    FramesSurroundLEDONMoving = FramesSurroundLEDONMoving + numSurroundON
                    VelocitySurroundMove = sum(VelocityColumn.loc[inSurroundIndexON])
                    TotalVelSurroundONMove = TotalVelSurroundONMove + VelocitySurroundMove
                    print("Frames Surround Move", FramesSurroundLEDONMoving)                

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
                print("LEDoffstart", currentLEDoffStart, 'LEDoffend', currentLEDoffEnd)
                MovementBlocksOFF +=1
                MovementDurationOFF = (currentLEDoffEnd-currentLEDoffStart)/30
                TotalMoveDurationOFF = TotalMoveDurationOFF + MovementDurationOFF
                TotalFramesMovingOFF = TotalFramesMovingOFF + numOFF
                MoveVelocityOFF = sum (VelocityColumn.loc[currentLEDoffStart:currentLEDoffEnd])
                TotalVelocityOFF = TotalVelocityOFF + MoveVelocityOFF
                isCenterSpanOFF = InZoneColumn.loc[currentLEDoffStart:currentLEDoffEnd]
                numCenterOFF = len(isCenterSpanOFF[isCenterSpanOFF == 1])
                numSurroundOFF = len(isCenterSpanOFF[isCenterSpanOFF == 0])
                inCenterIndexOFF = isCenterSpanOFF[isCenterSpanOFF == 1].index.tolist()
                inSurroundIndexOFF = isCenterSpanOFF[isCenterSpanOFF == 0].index.tolist()
                
                CenterThreshold = 0.98
                SurroundThreshold = 0.02
                PercentCenterOFF = numCenterOFF/numOFF
        
                if PercentCenterOFF >=CenterThreshold:
                    CenterLEDOFFOnly +=1
                    FramesCenterLEDOFFOnly = FramesCenterLEDOFFOnly + numCenterOFF
                    VelocityCenterOnly = sum(VelocityColumn.loc[inCenterIndexOFF])
                    TotalVelCenterOFFOnly = TotalVelCenterOFFOnly + VelocityCenterOnly
                    print("Frames Center Only", FramesCenterLEDOFFOnly)
                if PercentCenterOFF <=SurroundThreshold:
                    SurroundLEDOFFOnly +=1
                    FramesSurroundLEDOFFOnly = FramesSurroundLEDOFFOnly + numSurroundOFF
                    VelocitySurroundOnly = sum(VelocityColumn.loc[inSurroundIndexOFF])
                    TotalVelSurroundOFFOnly = TotalVelSurroundOFFOnly + VelocitySurroundOnly
                    print("Frames Surround Only", FramesSurroundLEDOFFOnly)
                if numCenterOFF > inCenterCutOff:
                    CenterLEDOFFTrue +=1
                    FramesCenterLEDOFFMoving = FramesCenterLEDOFFMoving + numCenterOFF
                    VelocityCenterMove = sum(VelocityColumn.loc[inCenterIndexOFF])
                    TotalVelCenterOFFMove = TotalVelCenterOFFMove + VelocityCenterMove
                    print("Frames Center Move", FramesCenterLEDOFFMoving)
                if numSurroundOFF > inCenterCutOff:
                    SurroundLEDOFFTrue +=1
                    FramesSurroundLEDOFFMoving = FramesSurroundLEDOFFMoving + numSurroundOFF
                    VelocitySurroundMove = sum(VelocityColumn.loc[inSurroundIndexOFF])
                    TotalVelSurroundOFFMove = TotalVelSurroundOFFMove + VelocitySurroundMove
                    print("Frames Surround Move", FramesSurroundLEDOFFMoving)
       
    PercentTimeinCenterMovingON = FramesCenterLEDONMoving/len(InZoneColumn)*100
    PercentTimeinCenterMovingOFF = FramesCenterLEDOFFMoving/len(InZoneColumn)*100
    PercentTimeinSurroundMovingON = FramesSurroundLEDONMoving/len(InZoneColumn)*100
    PercentTimeinSurroundMovingOFF = FramesSurroundLEDOFFMoving/len(InZoneColumn)*100
    PercentTimeCenterONOnly = FramesCenterLEDONOnly/len(InZoneColumn)*100
    PercentTimeSurroundONOnly = FramesSurroundLEDONOnly/len(InZoneColumn)*100
    PercentTimeCenterOFFOnly = FramesCenterLEDOFFOnly/len(InZoneColumn)*100
    PercentTimeSurroundOFFOnly = FramesSurroundLEDOFFOnly/len(InZoneColumn)*100
    AvgMoveDurationON = TotalMoveDurationON/MovementBlocksON
    AvgMoveDurationOFF = TotalMoveDurationOFF/MovementBlocksOFF
    AvgVelocityMovingON = TotalVelocityON/TotalFramesMovingON
    AvgVelocityMovingOFF = TotalVelocityOFF/TotalFramesMovingOFF
    if FramesCenterLEDONOnly > 0:
        AvgVelocityCenterONOnly = TotalVelCenterONOnly/FramesCenterLEDONOnly
    else:
        AvgVelocityCenterONOnly = 0
    if FramesSurroundLEDONOnly > 0:
        AvgVelocitySurroundONOnly = TotalVelSurroundONOnly/FramesSurroundLEDONOnly
    else:
        AvgVelocitySurroundONOnly = 0
    if FramesCenterLEDONMoving > 0:
        AvgVelocityCenterONMove = TotalVelCenterONMove/FramesCenterLEDONMoving
    else:
        AvgVelocityCenterONMove = 0
    if FramesSurroundLEDONMoving > 0:
        AvgVelocitySurroundONMove = TotalVelSurroundONMove/FramesSurroundLEDONMoving 
    else:
        AvgVelocitySurroundONMove = 0
    if FramesCenterLEDOFFOnly > 0:
        AvgVelocityCenterOFFOnly = TotalVelCenterOFFOnly/FramesCenterLEDOFFOnly
    else:
        AvgVelocityCenterOFFOnly = 0
    if FramesSurroundLEDOFFOnly > 0:
        AvgVelocitySurroundOFFOnly = TotalVelSurroundOFFOnly/FramesSurroundLEDOFFOnly
    else:
        AvgVelocitySurroundOFFOnly = 0
    if FramesCenterLEDOFFMoving > 0:
        AvgVelocityCenterOFFMove = TotalVelCenterOFFMove/FramesCenterLEDOFFMoving
    else:
        AvgVelocityCenterOFFMove = 0
    if FramesSurroundLEDOFFMoving > 0:
        AvgVelocitySurroundOFFMove = TotalVelSurroundOFFMove/FramesSurroundLEDOFFMoving 
    else:
        AvgVelocitySurroundOFFMove = 0

    return PercentTimeinCenterLEDon, PercentTimeinCenterLEDoff, PercentTimeinCenterMovingON, \
    PercentTimeinCenterMovingOFF, PercentTimeinSurroundMovingON, PercentTimeinSurroundMovingOFF,\
    PercentTimeCenterONOnly, PercentTimeCenterOFFOnly, PercentTimeSurroundONOnly, PercentTimeSurroundOFFOnly,\
    MovementBlocksON, MovementBlocksOFF, AvgMoveDurationON, AvgMoveDurationOFF, TotalMoveDurationON,\
    TotalMoveDurationOFF, AvgVelocityMovingON, AvgVelocityMovingOFF, AvgVelocityCenterONOnly,\
    AvgVelocityCenterOFFOnly, AvgVelocitySurroundONOnly, AvgVelocitySurroundOFFOnly, AvgVelocityCenterONMove,\
    AvgVelocityCenterOFFMove, AvgVelocitySurroundONMove, AvgVelocitySurroundOFFMove, \
    AvgVelocityLEDon, AvgVelocityLEDoff, CenterLEDONTrue, CenterLEDOFFTrue, SurroundLEDONTrue, \
    SurroundLEDOFFTrue, CenterLEDONOnly, CenterLEDOFFOnly, SurroundLEDONOnly, SurroundLEDOFFOnly
    
for File in FileList:    
    if not File.endswith('.xlsx'):
        print ("skipping file named", File)
        continue
    
    df = pd.read_excel(File, header = None)
    NumHead,Sbj,GT,LEDstate,DateTime,Note = Extract_Header_Info(df)
    NB10MIN_NHR = NumHead +1
    FiveMinBaseline2HSession_NHR = NumHead + 8992
    DataBlock = create_DataBlock(df, NO_BASELINE_10MIN_SESSION_TRIAL_LENGTH,NB10MIN_NHR,FiveMinBaseline2HSession_NHR)
    isMovingData = DataBlock.loc[:,'Movement(Moving / Center-point)']
    LEDon = DataBlock.loc[:,'LED ON']
    isInCenter = DataBlock.loc[:,'In zone']
    Velocity = DataBlock.loc [:, 'Velocity']
    startingPoints,endingPoints = Binary_Data_Transition_Point_Finder(isMovingData)
    if DataBlock['Trial time'].iloc[-1] > NO_BASELINE_10MIN_SESSION_TRIAL_LENGTH:
        startingPoints,endingPoints = Remove_Minute2_LED_OFF(DataBlock,startingPoints,endingPoints) 
        print ("2 minute LED off period file")
    PerTimeCenterON, PerTimeCenterOFF, PerTimeCenterMoveON,PerTimeCenterMoveOFF, PerTimeSurroundMoveON,\
    PerTimeSurroundMoveOFF, PerTimeCenterOnlyON,PerTimeCenterOnlyOFF, PerTimeSurroundOnlyON ,\
    PerTimeSurroundOnlyOFF, MovesON,MovesOFF, AvgMoveTimeON,AvgMoveTimeOFF, TotalMoveTimeON,TotalMoveTimeOFF,\
    AvgVelMoveON,AvgVelMoveOFF, AvgVelCO_ON, AvgVelCO_OFF, AvgVelSO_ON, AvgVelSO_OFF, AvgVelCM_ON, AvgVelCM_OFF,\
    AvgVelSM_ON, AvgVelSM_OFF,AvgVelON,AvgVelOFF,inCenterON, inCenterOFF,inSurroundON, inSurroundOFF,\
    inCenterOnlyON, inCenterOnlyOFF, InSurroundOnlyON, InSurroundOnlyOFF \
    = Movement_Analysis(DataBlock,isInCenter,Velocity,startingPoints, endingPoints)
   
    
    Dataz = [Sbj,GT, LEDstate, DateTime, Note, PerTimeCenterON, PerTimeCenterOFF, PerTimeCenterMoveON,\
    PerTimeCenterMoveOFF, PerTimeSurroundMoveON,PerTimeSurroundMoveOFF, PerTimeCenterOnlyON,\
    PerTimeCenterOnlyOFF, PerTimeSurroundOnlyON ,PerTimeSurroundOnlyOFF, MovesON,MovesOFF, AvgMoveTimeON,\
    AvgMoveTimeOFF, TotalMoveTimeON,TotalMoveTimeOFF,AvgVelMoveON,AvgVelMoveOFF, AvgVelCO_ON, \
    AvgVelCO_OFF, AvgVelSO_ON, AvgVelSO_OFF, AvgVelCM_ON, AvgVelCM_OFF,AvgVelSM_ON, AvgVelSM_OFF,\
    AvgVelON,AvgVelOFF,inCenterON, inCenterOFF,inSurroundON, inSurroundOFF,inCenterOnlyON, \
    inCenterOnlyOFF, InSurroundOnlyON, InSurroundOnlyOFF]
    myDataList.append(Dataz)
    print (myDataList)
    
resultsDF = pd.DataFrame(data = myDataList, columns=resultsColumns)
os.chdir("/Users/leblanckh/data/Opto_OpenField_RawData/OutputFiles")
writer = pd.ExcelWriter('Opto_OpenField_Movement_Analysis2.xlsx')
resultsDF.to_excel(writer,'Sheet1')
writer.save()





