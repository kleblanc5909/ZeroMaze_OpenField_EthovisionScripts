# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 21:28:53 2017

@author: leblanckh
"""

import pandas as pd
import numpy as np
import os
import sys
dataFolder = "/Users/leblanckh/data/KO_WT_ZeroMaze_RawData"
resultsColumns = ["Subject", "Group", "Gender", "Timestamp", "Notes", "Time in Open (%)", \
"Time in Open while moving (%)", "Time in Closed while moving (%)", "Time moving in Open only (%)",\
 "Time moving in Closed only (%)","Number of movements", "Average duration of movements (s)", \
 "Total time moving (s)", "Speed while moving (cm/s)", "Average Velocity for Open Only Movements (cm/s)",\
 "Average Velocity for Closed Only Movements (cm/s)", "Average Velocity for Open Movements (cm/s)", \
 "Average Velocity for Closed Movements (cm/s)", "Average Velocity (cm/s)", "Number of movements in Open",\
 "Number of movements in Closed", "Number of Open-only movements", "Number of Closed-only movements"]
myDataList = []
os.chdir(dataFolder)
FileList = os.listdir(dataFolder)

NO_BASELINE_10MIN_SESSION_TRIAL_LENGTH = 601
ONE_MINUTE_AS_FRAMES = 1798

def Extract_Header_Info(dataframe):
    """
    Identifies the length of the header,Extracts independent variable 
    information from the header, and assigns the column names to the row that
    contains the column information in the header
    
    inputs: dataframe
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
    
def Create_DataBlock(dataframe, cutoffTimeSeconds, shortTimeIdx, longTimeIdx):
    """
    Creates a dataFrame, theData, that contains the actual raw data values.
    It also replaces missing values in the movement column of data with NaN, then interpolates
    returns: theData, movement column data (isMovingData), LED on column (LEDon)
    """
    finalColumnName = dataframe.columns[-1]
    initialColumnName = dataframe.columns[0]

    if dataframe['Trial time'].iloc[-1] <cutoffTimeSeconds:
        theData = dataframe.loc[shortTimeIdx:,initialColumnName:finalColumnName] 
    else:
        theData = dataframe.loc[longTimeIdx:,initialColumnName:finalColumnName] 
    
    theData[:1].replace('-',0,inplace = True)
    theData[-1:].replace('-',0,inplace = True)
     #replace missing data (-) with NaN, then interpolate
    theData.replace('-',np.NaN,inplace = True)
    theData.interpolate(method = 'values', axis = 0, inplace = True)
    theData.loc[:,'Movement(Moving / center-point)'] = Binary_Data_Interpolation_CleanUp(theData.loc[:,'Movement(Moving / center-point)']) 
    theData.loc[:,'In open'] = Binary_Data_Interpolation_CleanUp(theData.loc[:,'In open']) 
    
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
    
def Movement_Analysis(InZoneColumn,VelocityColumn,MoveStart,MoveEnd):
    """
    Sets up the initial conditions, loops over the movement data and counts movement blocks,
    duration and velocity,identifies movements into the Open, and runs all of the calculations
    
    inputs: inZoneColumn, Velocity Column, and the starting and ending points of the moving column
    returns: PercentTimeinOpen, PercentTimeinOpenMoving,MovementBlocks, AvgMoveDuration, 
    TotalMoveDuration, AvgVelocityMoving, AvgVelocity,OpenTrue
    """
    OpenTrue = 0
    OpenOnly = 0
    ClosedOnly = 0
    ClosedTrue = 0
    inOpenCutOff = 15
    FramesOpenMoving = 0
    FramesOpenOnly = 0
    FramesClosedOnly = 0
    FramesClosedMoving = 0
    MovementBlocks = 0
    TotalMoveDuration = 0
    TotalVelocity = 0
    TotalVelOpenOnly = 0
    TotalVelClosedOnly = 0
    TotalVelOpenMove = 0
    TotalVelClosedMove = 0
    TotalFramesMoving = 0
    
    #loop over the movement data and identify if mouse enters Open
    if len(MoveStart) != len(MoveEnd):
        print ("Uneven start and end pairs. There are", (len(MoveStart)), "starting points and" , (len(MoveEnd)) ,"MoveEnd" )
    for i in range(len(MoveStart)):
        currentStart = MoveStart[i]
        currentEnd = MoveEnd[i] - 1
        #NOTE: MoveEnd are one frame beyond the end of the movement.  To correct for this, I have adjusted the current end point back one frame.
        MovementBlocks += 1
        MovementDuration = (currentEnd-currentStart)/30
        TotalMoveDuration = TotalMoveDuration + MovementDuration
        FramesMoving = (currentEnd-currentStart)
        TotalFramesMoving = TotalFramesMoving + FramesMoving
        MoveVelocity = sum (VelocityColumn.loc[currentStart:currentEnd])
        TotalVelocity = TotalVelocity + MoveVelocity
        isOpenSpan = InZoneColumn[currentStart:currentEnd]
        numOpen = len(isOpenSpan[isOpenSpan == 1])
        numClosed = len(isOpenSpan[isOpenSpan == 0])
        inOpenIndex = isOpenSpan[isOpenSpan == 1].index.tolist()
        inClosedIndex = isOpenSpan[isOpenSpan == 0].index.tolist()
        
        OpenThreshold = 0.95
        ClosedThreshold = 0.05
        PercentOpen = numOpen/FramesMoving
        
        if PercentOpen >=OpenThreshold:
            OpenOnly +=1
            FramesOpenOnly = FramesOpenOnly + numOpen
            VelocityOpenOnly = sum(VelocityColumn.loc[inOpenIndex])
            TotalVelOpenOnly = TotalVelOpenOnly + VelocityOpenOnly
        if PercentOpen <=ClosedThreshold:
            ClosedOnly +=1
            FramesClosedOnly = FramesClosedOnly + numClosed
            VelocityClosedOnly = sum(VelocityColumn.loc[inClosedIndex])
            TotalVelClosedOnly = TotalVelClosedOnly + VelocityClosedOnly
        if numOpen > inOpenCutOff:
            OpenTrue +=1
            FramesOpenMoving = FramesOpenMoving + numOpen
            VelocityOpenMove = sum(VelocityColumn.loc[inOpenIndex])
            TotalVelOpenMove = TotalVelOpenMove + VelocityOpenMove
        if numClosed > inOpenCutOff:
            ClosedTrue +=1
            FramesClosedMoving = FramesClosedMoving + numClosed
            VelocityClosedMove = sum(VelocityColumn.loc[inClosedIndex])
            TotalVelClosedMove = TotalVelClosedMove + VelocityClosedMove
            
    PercentTimeinOpenMoving = FramesOpenMoving/len(InZoneColumn)*100
    PercentTimeinClosedMoving = FramesClosedMoving/len(InZoneColumn)*100
    PercentTimeOpenOnly = FramesOpenOnly/len(InZoneColumn)*100
    PercentTimeClosedOnly = FramesClosedOnly/len(InZoneColumn)*100
    AvgMoveDuration = TotalMoveDuration/MovementBlocks
    AvgVelocityMoving = TotalVelocity/TotalFramesMoving
    if FramesOpenOnly > 0:
        AvgVelocityOpenOnly = TotalVelOpenOnly/FramesOpenOnly
    else:
        AvgVelocityOpenOnly = 0
    if FramesClosedOnly > 0:
        AvgVelocityClosedOnly = TotalVelClosedOnly/FramesClosedOnly
    else:
        AvgVelocityClosedOnly = 0
    if FramesOpenMoving > 0:
        AvgVelocityOpenMove = TotalVelOpenMove/FramesOpenMoving
    else:
        AvgVelocityOpenMove = 0
    if FramesClosedMoving > 0:
        AvgVelocityClosedMove = TotalVelClosedMove/FramesClosedMoving 
    else:
        AvgVelocityClosedMove = 0
    PercentTimeinOpen = sum (InZoneColumn[InZoneColumn == 1])/len(InZoneColumn)*100
    AvgVelocity = sum (VelocityColumn)/len(VelocityColumn)
    return PercentTimeinOpen, PercentTimeinOpenMoving, \
    PercentTimeinClosedMoving, PercentTimeOpenOnly, PercentTimeClosedOnly, \
    MovementBlocks, AvgMoveDuration,TotalMoveDuration, AvgVelocityMoving, \
    AvgVelocityOpenOnly, AvgVelocityClosedOnly, AvgVelocityOpenMove, AvgVelocityClosedMove,\
    AvgVelocity, OpenTrue, ClosedTrue, OpenOnly, ClosedOnly

for File in FileList:
    if not File.endswith('.xlsx'):
        print ("skipping file named", File)
        continue
    df = pd.read_excel(File, header = None)
    
    NumHead,Sbj,Gt,Sex,DateTime,Note = Extract_Header_Info(df)
    NB_10MIN = NumHead + 1
    ONE_MIN_BASE = NumHead + 1801
    DataBlock = Create_DataBlock(df, NO_BASELINE_10MIN_SESSION_TRIAL_LENGTH,NB_10MIN,ONE_MIN_BASE)
    isMovingData = DataBlock.loc[:,'Movement(Moving / center-point)']
    isInOpen = DataBlock.loc[:,'In open']
    Velocity = DataBlock.loc [:, 'Velocity']
    startingPoints,endingPoints = Binary_Data_Transition_Point_Finder(isMovingData)
    PerTimeOpen, PerTimeOpenMove,PerTimeClosedMove,PerTimeOpenOnly,PerTimeClosedOnly,\
    Moves,AvgMoveTime,TotalMoveTime,AvgVelMove,AvgVelCO, AvgVelSO, AvgVelCM, AvgVelSM,AvgVel,\
    inOpen, inClosed, inOpenOnly, InClosedOnly = Movement_Analysis(isInOpen,Velocity,startingPoints,endingPoints)
   
    Dataz = [Sbj,Gt, Sex, DateTime, Note, PerTimeOpen, PerTimeOpenMove,PerTimeClosedMove, \
    PerTimeOpenOnly, PerTimeClosedOnly, Moves, AvgMoveTime, TotalMoveTime, AvgVelMove, \
    AvgVelCO, AvgVelSO, AvgVelCM, AvgVelSM, AvgVel, inOpen, inClosed, inOpenOnly, InClosedOnly]
    myDataList.append(Dataz)
    print (myDataList)
  
resultsDF = pd.DataFrame(data = myDataList, columns=resultsColumns)
os.chdir("/Users/leblanckh/data/KO_WT_ZeroMaze_RawData/OutputFiles")
writer = pd.ExcelWriter('KO_and_WT_Zero_Maze_Movement_Analysis4.xlsx')
resultsDF.to_excel(writer,'Sheet1')
writer.save()


