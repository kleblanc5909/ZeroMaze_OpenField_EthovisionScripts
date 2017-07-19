# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 21:28:53 2017

@author: leblanckh
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import sys
import matplotlib.pyplot as plt
plt.style.use('ggplot')
dataFolder = "/Users/leblanckh/data/KO_WT_OpenField_RawData"
resultsColumns = ["Subject", "Group", "Gender", "Timestamp", "Notes", "Time in center (%)", "Time in center while moving (%)", "Number of movements", "Average duration of movements (s)", "Total time moving (s)", "Speed while moving (cm/s)", "Average Velocity (cm/s)", "Number of movements in center"]
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
    NumberofHeaderRows = int(dataframe.iloc[0,1])
    Header = dataframe.iloc[31:NumberofHeaderRows - 3,:]
    SubjectFilters = ["Subject", "Mouse","subject", "mouse"]
    Subject =  Header[Header[0].isin(SubjectFilters)].iloc[0,1]
    GenotypeFilters = ["Genotype", "Group", "genotype", "group", "genotype/group"]
    Genotype =  str(Header[Header[0].isin(GenotypeFilters)].iloc[0,1])
    GenderFilters = ["Sex", "Gender", "sex", "gender"]
    Gender =  Header[Header[0].isin(GenderFilters)].iloc[0,1]
    TimestampFilters = ["Timestamp", "timestamp"]
    Timestamp =  Header[Header[0].isin(TimestampFilters)].iloc[0,1]
    NotesFilters = ["Notes", "notes"]
    Notes =  str(Header[Header[0].isin(NotesFilters)].iloc[0,1])
    #print ("Subject:",Subject, "Genotype:", Genotype,"Gender:", Gender,"Timestamp:", Timestamp, "Notes", Notes)
    ColumnNames = dataframe.iloc[[NumberofHeaderRows - 2],:].values.tolist()
    dataframe.columns = ColumnNames
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
    
    firstRow = theData[:1]
    lastRow = theData[-1:]
    
    firstRow[firstRow == '-'] = 0
    lastRow[lastRow == '-'] = 0
    
    theData[:1] = firstRow
    theData[-1:] = lastRow
     #replace missing data (-) with NaN, then interpolate
    theData.replace('-',np.NaN,inplace = True)
    theData.interpolate(method = 'values', axis = 0, inplace = True)
    theData.loc[:,'Movement(Moving / Center-point)'] = Binary_Data_Interpolation_CleanUp(theData.loc[:,'Movement(Moving / Center-point)']) 
    theData.loc[:,'In zone'] = Binary_Data_Interpolation_CleanUp(theData.loc[:,'In zone']) 
    
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
    duration and velocity,identifies movements into the center, and runs all of the calculations
    
    inputs: inZoneColumn, Velocity Column, and the starting and ending points of the moving column
    returns: PercentTimeinCenter, PercentTimeinCenterMoving,MovementBlocks, AvgMoveDuration, 
    TotalMoveDuration, AvgVelocityMoving, AvgVelocity,CenterTrue
    """
    CenterTrue = 0
    inCenterCutOff = 15
    FramesCenterMoving = 0
    MovementBlocks = 0
    TotalMoveDuration = 0
    TotalVelocity = 0
    TotalFramesMoving = 0

        
    #loop over the movement data and identify if mouse enters Center
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
    PercentTimeinCenter = sum (isInCenter[isInCenter == 1])/len(isInCenter)*100
    AvgVelocity = sum (Velocity)/len(Velocity)
    return PercentTimeinCenter, PercentTimeinCenterMoving, MovementBlocks, AvgMoveDuration, \
    TotalMoveDuration, AvgVelocityMoving, AvgVelocity, CenterTrue

def Create_Group_Bar_Plot (dataframe,dataColumnName):
    """
    Groups the dataframe by Genotype and one of the measurements, 
    takes the mean and standard deviation, and plots the data on a bar graph
    
    input: dataframe, dataColumnName
    return:
    """
    GroupedResults = dataframe.groupby('Group')[dataColumnName]
    print("GroupedResults", GroupedResults)
    GroupMeans = GroupedResults.mean()
    print("groupMeans", GroupMeans)
    GroupError = GroupedResults.std()
    #color_list = ['tab:gray', 'tab:red', 'tab:gray', 'tab:orange', 'tab:gray', 'tab:purple']
    fig,ax = plt.subplots()
    ax.set_ylabel(dataColumnName)
    GroupMeans.plot.bar(yerr=GroupError, ax = ax)
    
for File in FileList:
    if not File.endswith('.xlsx'):
        print ("skipping file named", File)
        continue
    df = pd.read_excel(File, header = None)
    
    NumHead,Sbj,Gt,Sex,DateTime,Note = Extract_Header_Info(df)
    NB_10MIN = NumHead + 1
    ONE_MIN_BASE = NumHead + 1801
    DataBlock = Create_DataBlock(df, NO_BASELINE_10MIN_SESSION_TRIAL_LENGTH,NB_10MIN,ONE_MIN_BASE)
    isMovingData = DataBlock.loc[:,'Movement(Moving / Center-point)']
    isInCenter = DataBlock.loc[:,'In zone']
    Velocity = DataBlock.loc [:, 'Velocity']
    startingPoints,endingPoints = Binary_Data_Transition_Point_Finder(isMovingData)
    PerTimeCenter, PerTimeCenterMove,Moves,AvgMoveTime,TotalMoveTime,AvgVelMove,AvgVel,inCenter \
    = Movement_Analysis(isInCenter,Velocity,startingPoints,endingPoints)
   
    Dataz = [Sbj,Gt, Sex, DateTime, Note, PerTimeCenter, PerTimeCenterMove, \
    Moves, AvgMoveTime, TotalMoveTime, AvgVelMove, AvgVel, \
    inCenter]
    myDataList.append(Dataz)
    print (myDataList)
    
resultsDF = pd.DataFrame(data = myDataList, columns=resultsColumns)
ColumnsToAnalyze = resultsColumns[5:]
for ColumnLabel in ColumnsToAnalyze:
    Create_Group_Bar_Plot(resultsDF,ColumnLabel)

os.chdir("/Users/leblanckh/data/KO_WT_OpenField_RawData/OutputFiles")
writer = pd.ExcelWriter('KO_and_WT_OpenField_Movement_Analysis.xlsx')
resultsDF.to_excel(writer,'Sheet1')
writer.save()


