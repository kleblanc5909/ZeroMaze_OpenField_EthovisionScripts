# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 21:28:53 2017

@author: leblanckh
"""

import pandas as pd
import numpy as np
import os
import sys
dataFolder = "/Users/leblanckh/data/Opto_ZeroMaze_RawData_Combined"
resultsColumns = ["Subject", "Group", "LED Power", "Timestamp", "Notes", "Time in Open, LED on (%)",\
"Time in Open, LED off (%)",  "Time in Open while moving and LED on(%)", "Time in Open while moving and LED off(%)",\
"Time in Closed while moving and LED on(%)", "Time in Closed while moving and LED off(%)", \
"Time moving in Open only, LED on (%)", "Time moving in Open only, LED off (%)",\
"Time moving in Closed only, LED on (%)", "Time moving in Closed only, LED off (%)",\
"Number of movements, LED on", "Number of movements, LED off","Average duration of movements, LED on (s)",\
 "Average duration of movements, LED off (s)", "Total time moving, LED on (s)", "Total time moving, LED off (s)",\
 "Speed while moving, LED on (cm/s)", "Speed while moving, LED off (cm/s)",\
 "Average Velocity for Open Only Movements, LED on (cm/s)","Average Velocity for Open Only Movements, LED off (cm/s)", \
 "Average Velocity for Closed Only Movements, LED on (cm/s)", "Average Velocity for Closed Only Movements, LED off (cm/s)", \
 "Average Velocity for Open Movements, LED on (cm/s)","Average Velocity for Open Movements, LED off (cm/s)", \
 "Average Velocity for Closed Movements, LED on (cm/s)","Average Velocity for Closed Movements, LED off (cm/s)", \
 "Average velocity, LED on (cm/s)","Average velocity, LED off (cm/s)", \
 "Number of movements in Open, LED on", "Number of movements in Open, LED off", \
 "Number of movements in Closed, LED on","Number of movements in Closed, LED off",\
 "Number of Open-only movements, LED on", "Number of Open-only movements, LED off", \
 "Number of Closed-only movements, LED on", "Number of Closed-only movements, LED off"]
myDataList = []
os.chdir(dataFolder)
FileList = os.listdir(dataFolder)

ONE_MINUTE_AS_FRAMES = 1798
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

def Create_DataBlock(dataframe, TimeIdx):
    """
    Creates a dataFrame, theData, that contains the actual raw data values.
    It also replaces missing values in the movement column of data with NaN, then interpolates
    returns: theData, movement column data (isMovingData), LED on column (LEDon)
    """
    finalColumnName = dataframe.columns[-1]
    initialColumnName = dataframe.columns[0]
    
    theData = dataframe.loc[TimeIdx:,initialColumnName:finalColumnName] #NumHdr + 8992
    
    theData[:1].replace('-',0,inplace = True)
    theData[-1:].replace('-',0,inplace = True)
     #replace missing data (-) with NaN, then interpolate
    theData.replace('-',np.NaN,inplace = True)
    theData.interpolate(method = 'values', axis = 0, inplace = True)
    theData.loc[:,'Movement(Moving / Open-point)'] = Binary_Data_Interpolation_CleanUp(theData.loc[:,'Movement(Moving / Center-point)']) 
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
    
def Movement_Analysis(theData, InZoneColumn,VelocityColumn,MoveStart,MoveEnd):
    """
    Sets up the initial conditions, loops over the movement data and counts movement blocks,
    duration and velocity,identifies movements into the Open, and runs all of the calculations
    
    
    inputs: inZoneColumn, Velocity Column, and the starting and ending points of the moving column
    returns: PercentTimeinOpen, PercentTimeinOpenMoving,MovementBlocks, AvgMoveDuration, 
    TotalMoveDuration, AvgVelocityMoving, AvgVelocity,OpenTrue
    """
    
    #isMovingData = theData.loc[:,'Movement(Moving / Center-point)']
    InZoneColumn = theData.loc[:,'In zone']
    LEDon = theData.loc[:,'LED ON']
    #TotalNumberDataRows = len(LEDon)
    LEDoff = theData.loc[:,'LED OFF']
    VelocityColumn = theData.loc [:, 'Velocity']
    #PercentTimeinOpen = sum (InZoneColumn[InZoneColumn == 1])/len(InZoneColumn)*100
    #AvgVelocity = sum (VelocityColumn)/len(VelocityColumn)
    
    LEDonIdx = LEDon[LEDon ==1].index.tolist()
    LEDon_Filtered_InOpen = InZoneColumn[LEDonIdx]
    LEDon_Filtered_Velocity = VelocityColumn[LEDonIdx]
    PercentTimeinOpenLEDon = sum (LEDon_Filtered_InOpen[0:])/len(LEDonIdx)*100
    AvgVelocityLEDon = sum(LEDon_Filtered_Velocity)/len(LEDonIdx)
    LEDoffIdx = LEDoff[LEDoff ==1].index.tolist()
    LEDoff_Filtered_InOpen = InZoneColumn[LEDoffIdx]
    LEDoff_Filtered_Velocity = VelocityColumn[LEDoffIdx]
    PercentTimeinOpenLEDoff = sum (LEDoff_Filtered_InOpen[0:])/len(LEDoffIdx)*100
    AvgVelocityLEDoff = sum(LEDoff_Filtered_Velocity)/len(LEDoffIdx)
    
    
    OpenLEDONTrue = 0
    OpenLEDOFFTrue = 0
    OpenLEDONOnly = 0
    ClosedLEDONOnly = 0
    ClosedLEDONTrue = 0
    OpenLEDOFFOnly = 0
    ClosedLEDOFFOnly = 0
    ClosedLEDOFFTrue = 0
    inOpenCutOff = 15
    LEDCutOff = 15
    FramesOpenLEDONMoving = 0
    FramesOpenLEDOFFMoving = 0
    FramesOpenLEDONOnly = 0
    FramesClosedLEDONOnly = 0
    FramesClosedLEDONMoving = 0
    FramesOpenLEDOFFOnly = 0
    FramesClosedLEDOFFOnly = 0
    FramesClosedLEDOFFMoving = 0
    MovementBlocksON = 0
    TotalMoveDurationON = 0
    TotalVelocityON = 0
    TotalVelOpenONOnly = 0
    TotalVelClosedONOnly = 0
    TotalVelOpenONMove = 0
    TotalVelClosedONMove = 0
    TotalVelOpenOFFOnly = 0
    TotalVelClosedOFFOnly = 0
    TotalVelOpenOFFMove = 0
    TotalVelClosedOFFMove = 0
    TotalFramesMovingON = 0
    MovementBlocksOFF = 0
    TotalMoveDurationOFF = 0
    TotalVelocityOFF = 0
    TotalFramesMovingOFF = 0

    
    #loop over the movement data and identify occurances of mouse entering Open
    if len(MoveStart) != len(MoveEnd):
        print ("Uneven start and end pairs. There are", (len(MoveStart)), "starting points and" , (len(MoveEnd)) ,"endingPoints" )
    for i in range(len(MoveStart)):
        currentStart = MoveStart[i]+1
        currentEnd = MoveEnd[i]
#        print("current start, end", currentStart, currentEnd)
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
#                print ("current ON start:", currentLEDonStart, "current ON end:", currentLEDonEnd)
                MovementBlocksON +=1
                MovementDurationON = (currentLEDonEnd-currentLEDonStart)/30
                TotalMoveDurationON = TotalMoveDurationON + MovementDurationON
                TotalFramesMovingON = TotalFramesMovingON + numON
                MoveVelocityON = sum (VelocityColumn.loc[currentLEDonStart:currentLEDonEnd])
                TotalVelocityON = TotalVelocityON + MoveVelocityON
                isOpenSpanON = InZoneColumn.loc[currentLEDonStart:currentLEDonEnd]
                numOpenON = len(isOpenSpanON[isOpenSpanON == 1])
                numClosedON = len(isOpenSpanON[isOpenSpanON == 0])
                inOpenIndexON = isOpenSpanON[isOpenSpanON == 1].index.tolist()
                inClosedIndexON = isOpenSpanON[isOpenSpanON == 0].index.tolist()
                
                OpenThreshold = 0.98
                ClosedThreshold = 0.02
                PercentOpenON = numOpenON/numON
        
                if PercentOpenON >=OpenThreshold:
                    OpenLEDONOnly +=1
                    FramesOpenLEDONOnly = FramesOpenLEDONOnly + numOpenON
                    VelocityOpenOnly = sum(VelocityColumn.loc[inOpenIndexON])
                    TotalVelOpenONOnly = TotalVelOpenONOnly + VelocityOpenOnly
#                    print("Frames Open Only", FramesOpenLEDONOnly)
                if PercentOpenON <=ClosedThreshold:
                    ClosedLEDONOnly +=1
                    FramesClosedLEDONOnly = FramesClosedLEDONOnly + numClosedON
                    VelocityClosedOnly = sum(VelocityColumn.loc[inClosedIndexON])
                    TotalVelClosedONOnly = TotalVelClosedONOnly + VelocityClosedOnly
#                    print("Frames Closed Only", FramesClosedLEDONOnly)
                if numOpenON > inOpenCutOff:
                    OpenLEDONTrue +=1
                    FramesOpenLEDONMoving = FramesOpenLEDONMoving + numOpenON
                    VelocityOpenMove = sum(VelocityColumn.loc[inOpenIndexON])
                    TotalVelOpenONMove = TotalVelOpenONMove + VelocityOpenMove
#                    print("Frames Open Move", FramesOpenLEDONMoving)
                if numClosedON > inOpenCutOff:
                    ClosedLEDONTrue +=1
                    FramesClosedLEDONMoving = FramesClosedLEDONMoving + numClosedON
                    VelocityClosedMove = sum(VelocityColumn.loc[inClosedIndexON])
                    TotalVelClosedONMove = TotalVelClosedONMove + VelocityClosedMove
#                    print("Frames Closed Move", FramesClosedLEDONMoving)                

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
#                print("LEDoffstart", currentLEDoffStart, 'LEDoffend', currentLEDoffEnd)
                MovementBlocksOFF +=1
                MovementDurationOFF = (currentLEDoffEnd-currentLEDoffStart)/30
                TotalMoveDurationOFF = TotalMoveDurationOFF + MovementDurationOFF
                TotalFramesMovingOFF = TotalFramesMovingOFF + numOFF
                MoveVelocityOFF = sum (VelocityColumn.loc[currentLEDoffStart:currentLEDoffEnd])
                TotalVelocityOFF = TotalVelocityOFF + MoveVelocityOFF
                isOpenSpanOFF = InZoneColumn.loc[currentLEDoffStart:currentLEDoffEnd]
                numOpenOFF = len(isOpenSpanOFF[isOpenSpanOFF == 1])
                numClosedOFF = len(isOpenSpanOFF[isOpenSpanOFF == 0])
                inOpenIndexOFF = isOpenSpanOFF[isOpenSpanOFF == 1].index.tolist()
                inClosedIndexOFF = isOpenSpanOFF[isOpenSpanOFF == 0].index.tolist()
                
                OpenThreshold = 0.98
                ClosedThreshold = 0.02
                PercentOpenOFF = numOpenOFF/numOFF
        
                if PercentOpenOFF >=OpenThreshold:
                    OpenLEDOFFOnly +=1
                    FramesOpenLEDOFFOnly = FramesOpenLEDOFFOnly + numOpenOFF
                    VelocityOpenOnly = sum(VelocityColumn.loc[inOpenIndexOFF])
                    TotalVelOpenOFFOnly = TotalVelOpenOFFOnly + VelocityOpenOnly
#                    print("Frames Open Only", FramesOpenLEDOFFOnly)
                if PercentOpenOFF <=ClosedThreshold:
                    ClosedLEDOFFOnly +=1
                    FramesClosedLEDOFFOnly = FramesClosedLEDOFFOnly + numClosedOFF
                    VelocityClosedOnly = sum(VelocityColumn.loc[inClosedIndexOFF])
                    TotalVelClosedOFFOnly = TotalVelClosedOFFOnly + VelocityClosedOnly
#                    print("Frames Closed Only", FramesClosedLEDOFFOnly)
                if numOpenOFF > inOpenCutOff:
                    OpenLEDOFFTrue +=1
                    FramesOpenLEDOFFMoving = FramesOpenLEDOFFMoving + numOpenOFF
                    VelocityOpenMove = sum(VelocityColumn.loc[inOpenIndexOFF])
                    TotalVelOpenOFFMove = TotalVelOpenOFFMove + VelocityOpenMove
#                    print("Frames Open Move", FramesOpenLEDOFFMoving)
                if numClosedOFF > inOpenCutOff:
                    ClosedLEDOFFTrue +=1
                    FramesClosedLEDOFFMoving = FramesClosedLEDOFFMoving + numClosedOFF
                    VelocityClosedMove = sum(VelocityColumn.loc[inClosedIndexOFF])
                    TotalVelClosedOFFMove = TotalVelClosedOFFMove + VelocityClosedMove
#                    print("Frames Closed Move", FramesClosedLEDOFFMoving)
       
    PercentTimeinOpenMovingON = FramesOpenLEDONMoving/len(InZoneColumn)*100
    PercentTimeinOpenMovingOFF = FramesOpenLEDOFFMoving/len(InZoneColumn)*100
    PercentTimeinClosedMovingON = FramesClosedLEDONMoving/len(InZoneColumn)*100
    PercentTimeinClosedMovingOFF = FramesClosedLEDOFFMoving/len(InZoneColumn)*100
    PercentTimeOpenONOnly = FramesOpenLEDONOnly/len(InZoneColumn)*100
    PercentTimeClosedONOnly = FramesClosedLEDONOnly/len(InZoneColumn)*100
    PercentTimeOpenOFFOnly = FramesOpenLEDOFFOnly/len(InZoneColumn)*100
    PercentTimeClosedOFFOnly = FramesClosedLEDOFFOnly/len(InZoneColumn)*100
    AvgMoveDurationON = TotalMoveDurationON/MovementBlocksON
    AvgMoveDurationOFF = TotalMoveDurationOFF/MovementBlocksOFF
    AvgVelocityMovingON = TotalVelocityON/TotalFramesMovingON
    AvgVelocityMovingOFF = TotalVelocityOFF/TotalFramesMovingOFF
    if FramesOpenLEDONOnly > 0:
        AvgVelocityOpenONOnly = TotalVelOpenONOnly/FramesOpenLEDONOnly
    else:
        AvgVelocityOpenONOnly = 0
    if FramesClosedLEDONOnly > 0:
        AvgVelocityClosedONOnly = TotalVelClosedONOnly/FramesClosedLEDONOnly
    else:
        AvgVelocityClosedONOnly = 0
    if FramesOpenLEDONMoving > 0:
        AvgVelocityOpenONMove = TotalVelOpenONMove/FramesOpenLEDONMoving
    else:
        AvgVelocityOpenONMove = 0
    if FramesClosedLEDONMoving > 0:
        AvgVelocityClosedONMove = TotalVelClosedONMove/FramesClosedLEDONMoving 
    else:
        AvgVelocityClosedONMove = 0
    if FramesOpenLEDOFFOnly > 0:
        AvgVelocityOpenOFFOnly = TotalVelOpenOFFOnly/FramesOpenLEDOFFOnly
    else:
        AvgVelocityOpenOFFOnly = 0
    if FramesClosedLEDOFFOnly > 0:
        AvgVelocityClosedOFFOnly = TotalVelClosedOFFOnly/FramesClosedLEDOFFOnly
    else:
        AvgVelocityClosedOFFOnly = 0
    if FramesOpenLEDOFFMoving > 0:
        AvgVelocityOpenOFFMove = TotalVelOpenOFFMove/FramesOpenLEDOFFMoving
    else:
        AvgVelocityOpenOFFMove = 0
    if FramesClosedLEDOFFMoving > 0:
        AvgVelocityClosedOFFMove = TotalVelClosedOFFMove/FramesClosedLEDOFFMoving 
    else:
        AvgVelocityClosedOFFMove = 0

    return PercentTimeinOpenLEDon, PercentTimeinOpenLEDoff, PercentTimeinOpenMovingON, \
    PercentTimeinOpenMovingOFF, PercentTimeinClosedMovingON, PercentTimeinClosedMovingOFF,\
    PercentTimeOpenONOnly, PercentTimeOpenOFFOnly, PercentTimeClosedONOnly, PercentTimeClosedOFFOnly,\
    MovementBlocksON, MovementBlocksOFF, AvgMoveDurationON, AvgMoveDurationOFF, TotalMoveDurationON,\
    TotalMoveDurationOFF, AvgVelocityMovingON, AvgVelocityMovingOFF, AvgVelocityOpenONOnly,\
    AvgVelocityOpenOFFOnly, AvgVelocityClosedONOnly, AvgVelocityClosedOFFOnly, AvgVelocityOpenONMove,\
    AvgVelocityOpenOFFMove, AvgVelocityClosedONMove, AvgVelocityClosedOFFMove, \
    AvgVelocityLEDon, AvgVelocityLEDoff, OpenLEDONTrue, OpenLEDOFFTrue, ClosedLEDONTrue, \
    ClosedLEDOFFTrue, OpenLEDONOnly, OpenLEDOFFOnly, ClosedLEDONOnly, ClosedLEDOFFOnly

for File in FileList:
    if not File.endswith('.xlsx'):
        print ("skipping file named", File)
        continue
    df = pd.read_excel(File, header = None)
    NumHead,Sbj,GT,LEDstate,DateTime,Note = Extract_Header_Info(df)
    ONE_MIN_BASE = NumHead + 1801
    DataBlock = Create_DataBlock(df,ONE_MIN_BASE)
    isMovingData = DataBlock.loc[:,'Movement(Moving / Open-point)']    
    isInOpen = DataBlock.loc[:,'In zone']
    LEDon = DataBlock.loc[:,'LED ON']
    LEDoff = DataBlock.loc[:,'LED OFF']
    Velocity = DataBlock.loc [:, 'Velocity']
    startingPoints,endingPoints = Binary_Data_Transition_Point_Finder(isMovingData)
    PerTimeOpenON, PerTimeOpenOFF, PerTimeOpenMoveON,PerTimeOpenMoveOFF, PerTimeClosedMoveON,\
    PerTimeClosedMoveOFF, PerTimeOpenOnlyON,PerTimeOpenOnlyOFF, PerTimeClosedOnlyON ,\
    PerTimeClosedOnlyOFF, MovesON,MovesOFF, AvgMoveTimeON,AvgMoveTimeOFF, TotalMoveTimeON,TotalMoveTimeOFF,\
    AvgVelMoveON,AvgVelMoveOFF, AvgVelCO_ON, AvgVelCO_OFF, AvgVelSO_ON, AvgVelSO_OFF, AvgVelCM_ON, AvgVelCM_OFF,\
    AvgVelSM_ON, AvgVelSM_OFF,AvgVelON,AvgVelOFF,inOpenON, inOpenOFF,inClosedON, inClosedOFF,\
    inOpenOnlyON, inOpenOnlyOFF, InClosedOnlyON, InClosedOnlyOFF \
    = Movement_Analysis(DataBlock,isInOpen,Velocity,startingPoints, endingPoints)

   
    Dataz = [Sbj,GT, LEDstate, DateTime, Note, PerTimeOpenON, PerTimeOpenOFF, PerTimeOpenMoveON,\
    PerTimeOpenMoveOFF, PerTimeClosedMoveON,PerTimeClosedMoveOFF, PerTimeOpenOnlyON,\
    PerTimeOpenOnlyOFF, PerTimeClosedOnlyON ,PerTimeClosedOnlyOFF, MovesON,MovesOFF, AvgMoveTimeON,\
    AvgMoveTimeOFF, TotalMoveTimeON,TotalMoveTimeOFF,AvgVelMoveON,AvgVelMoveOFF, AvgVelCO_ON, \
    AvgVelCO_OFF, AvgVelSO_ON, AvgVelSO_OFF, AvgVelCM_ON, AvgVelCM_OFF,AvgVelSM_ON, AvgVelSM_OFF,\
    AvgVelON,AvgVelOFF,inOpenON, inOpenOFF,inClosedON, inClosedOFF,inOpenOnlyON, \
    inOpenOnlyOFF, InClosedOnlyON, InClosedOnlyOFF]
    myDataList.append(Dataz)
    print (myDataList)
    
resultsDF = pd.DataFrame(data = myDataList, columns=resultsColumns)
os.chdir("/Users/leblanckh/data/Opto_ZeroMaze_RawData/OutputFiles")
writer = pd.ExcelWriter('Opto_Zero_Maze_Movement_Analysis3.xlsx')
resultsDF.to_excel(writer,'Sheet1')
writer.save()


