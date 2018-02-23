#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 10:34:49 2018

@author: thomas
"""

import numpy as np
import numpy.ma as ma
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import os
import re

#I would like to try to make this script such that ALL the information is 
# provided and I do not need to change parameters

#Make it more general! make it so that it can read in all the files from all the directories
#1) Loop over all directories in list
#2) Obtain hier_data files from each directory
#3) Read files and calulcate net force
#4) Create separate plot of net force for each case (like before) but now with COM!
#5) End of loop: Compile all net forces for each directory in array
#6) Plot force vectors for 1) Upper 2) Lower 3) COM

#Calculate the net force on our dimer through an oscillation
#1) Calculate the net force on each vertex point
#2) Average the force per timestep (weighted equally?) and use to determine 
#   net force on dimer through oscillation

PERIOD = 0.125

def getFilesForces(dirName,nvert,nDimers,Re,comArr,forceArr):
    #Data Dump Number
    timeFinal = 5.0
    dt = 1.0e-05
    dumpInterval = 100
    maxFileNumber = int(timeFinal/(dumpInterval*dt))*100
    fileNumber = np.arange(100,maxFileNumber+dumpInterval,dumpInterval)
    print('fileNumber[10] = ',fileNumber[10])
    
    #Allocate Arrays
    #Force on vertices
    #1)
    forceBotUp1Arr = np.zeros((2,nvert[0,0],len(fileNumber)))
    forceBotUp2Arr = np.zeros((2,nvert[0,1],len(fileNumber)))
    forceBotLow1Arr = np.zeros((2,nvert[1,0],len(fileNumber)))
    forceBotLow2Arr = np.zeros((2,nvert[1,1],len(fileNumber)))
    #Net Force per timestep
    #2)
    netForceBotUp1 = np.zeros((2,len(fileNumber)))
    netForceBotUp2 = np.zeros((2,len(fileNumber)))
    netForceBotLow1 = np.zeros((2,len(fileNumber)))
    netForceBotLow2 = np.zeros((2,len(fileNumber)))
    #Net Force During a Period
    netForcePeriodBotUp1 = np.zeros((2,int(timeFinal/PERIOD)))
    netForcePeriodBotUp2 = np.zeros((2,int(timeFinal/PERIOD)))
    netForcePeriodBotLow1 = np.zeros((2,int(timeFinal/PERIOD)))
    netForcePeriodBotLow2 = np.zeros((2,int(timeFinal/PERIOD)))
    
    for fN in fileNumber:
        fNindx = int(fN/dumpInterval)-1
        if(dirName == '../10/'):
            print(fNindx+1)
        if(fN < 1000):
            forceFilePath = dirName+'F.00'+str(fN)
        elif(fN >= 1000 and fN < 10000):
            forceFilePath = dirName+'F.0'+str(fN)
        else:
            forceFilePath = dirName+'F.'+str(fN)
        #Read force file line by line and save the values in a list using the \n delimiter
        lines = [line.strip() for line in open(forceFilePath)]
        #Break each list element into an array of numbers using the space delimiter
        lines = [line.split() for line in lines]
        
        lineNumber = 0
        vertNumber = 0
        while(vertNumber < nvert[0,0] or lineNumber > 2*nvert[0,0]+20):
            if(re.search('[a-zA-Z]',str(lines[lineNumber])) and 'e-0' not in str(lines[lineNumber])):
            #if (re.search('[a-zA-Z]',str(lines[lineNumber])) ):
                #print('Not a Number! line = ',lines[lineNumber])
                lineNumber+=1
            else:
                #print(float(lines[lineNumber][0]))
                #print(float(lines[lineNumber+1][0]))
                forceBotUp1Arr[0,vertNumber,fNindx] = float(lines[lineNumber][0])
                forceBotUp1Arr[1,vertNumber,fNindx] = float(lines[lineNumber+1][0])
                #print(forceBotUp1Arr[:,vertNumber,fN])
                lineNumber+=2
                vertNumber+=1
        #print('End of BotUp1')
        vertNumber = 0
        while(vertNumber < nvert[1,0] or lineNumber > 2*(nvert[0,0] + nvert[1,0])+20):
            if (re.search('[a-zA-Z]',str(lines[lineNumber])) and 'e-0' not in str(lines[lineNumber])):
                #print('Not a Number! line = ',lines[lineNumber])
                lineNumber+=1
            else:
                #print(float(lines[lineNumber][0]))
                #print(float(lines[lineNumber+1][0]))
                forceBotLow1Arr[0,vertNumber,fNindx] = float(lines[lineNumber][0])
                forceBotLow1Arr[1,vertNumber,fNindx] = float(lines[lineNumber+1][0])
                #print(forceBotUp1Arr[:,vertNumber,fN])
                lineNumber+=2
                vertNumber+=1
        vertNumber = 0
        #print('End of BotLow1')
        while(vertNumber < nvert[0,1] or lineNumber > 2*(nvert[0,0] + nvert[0,1]+nvert[1,0])+20):
            if (re.search('[a-zA-Z]',str(lines[lineNumber])) and 'e-0' not in str(lines[lineNumber])):
                #print('Not a Number! line = ',lines[lineNumber])
                lineNumber+=1
            else:
                #print(float(lines[lineNumber][0]))
                #print(float(lines[lineNumber+1][0]))
                forceBotUp2Arr[0,vertNumber,fNindx] = float(lines[lineNumber][0])
                forceBotUp2Arr[1,vertNumber,fNindx] = float(lines[lineNumber+1][0])
                #print(forceBotUp1Arr[:,vertNumber,fN])
                lineNumber+=2
                vertNumber+=1
        vertNumber = 0
        #print('End of BotUp2')
        while(vertNumber < nvert[1,1] or lineNumber > 2*(np.sum(nvert))+20):
            if (re.search('[a-zA-Z]',str(lines[lineNumber])) and 'e-0' not in str(lines[lineNumber])):
                #print('Not a Number! line = ',lines[lineNumber])
                lineNumber+=1
            else:
                #print(float(lines[lineNumber][0]))
                #print(float(lines[lineNumber+1][0]))
                forceBotLow2Arr[0,vertNumber,fNindx] = float(lines[lineNumber][0])
                forceBotLow2Arr[1,vertNumber,fNindx] = float(lines[lineNumber+1][0])
                #print(forceBotUp1Arr[:,vertNumber,fN])
                lineNumber+=2
                vertNumber+=1
        #print('End of BotLow2')
        #print(lines[lineNumber])
        #print(lines[lineNumber+1])
        vertNumber = 0
        
        #2) Calculate the Average Force for that TimeStep (All 4 structures)
        netForceBotUp1[0,fNindx] = np.sum(forceBotUp1Arr[0,:,fNindx])/float(nvert[0,0])
        netForceBotUp1[1,fNindx] = np.sum(forceBotUp1Arr[1,:,fNindx])/float(nvert[0,0])
        netForceBotUp2[0,fNindx] = np.sum(forceBotUp2Arr[0,:,fNindx])/float(nvert[0,1])
        netForceBotUp2[1,fNindx] = np.sum(forceBotUp2Arr[1,:,fNindx])/float(nvert[0,1])
        netForceBotLow1[0,fNindx] = np.sum(forceBotLow1Arr[0,:,fNindx])/float(nvert[1,0])
        netForceBotLow1[1,fNindx] = np.sum(forceBotLow1Arr[1,:,fNindx])/float(nvert[1,0])
        netForceBotLow2[0,fNindx] = np.sum(forceBotLow2Arr[0,:,fNindx])/float(nvert[1,1])
        netForceBotLow2[1,fNindx] = np.sum(forceBotLow2Arr[1,:,fNindx])/float(nvert[1,1])

    for i in range(int(timeFinal/PERIOD)):
        #Calculate Net Force on Dimer Piece!
        netForcePeriodBotUp1[0,i] = np.sum(netForceBotUp1[0,125*i:125*(i+1)])/(PERIOD/(dumpInterval*dt))
        netForcePeriodBotUp2[0,i] = np.sum(netForceBotUp2[0,125*i:125*(i+1)])/(PERIOD/(dumpInterval*dt))
        netForcePeriodBotLow1[0,i] = np.sum(netForceBotLow1[0,125*i:125*(i+1)])/(PERIOD/(dumpInterval*dt))
        netForcePeriodBotLow2[0,i] = np.sum(netForceBotLow2[0,125*i:125*(i+1)])/(PERIOD/(dumpInterval*dt))
        netForcePeriodBotUp1[1,i] = np.sum(netForceBotUp1[1,125*i:125*(i+1)])/(PERIOD/(dumpInterval*dt))
        netForcePeriodBotUp2[1,i] = np.sum(netForceBotUp2[1,125*i:125*(i+1)])/(PERIOD/(dumpInterval*dt))
        netForcePeriodBotLow1[1,i] = np.sum(netForceBotLow1[1,125*i:125*(i+1)])/(PERIOD/(dumpInterval*dt))
        netForcePeriodBotLow2[1,i] = np.sum(netForceBotLow2[1,125*i:125*(i+1)])/(PERIOD/(dumpInterval*dt))

    #Plot Net Force per Oscillation
    timeOsc = np.arange(0.0,timeFinal,PERIOD)
    plt.figure(0)
    plt.title('Net Force BotUp1: Re'+str(Re))
    plt.plot(timeOsc,netForcePeriodBotUp1[1,:],label = 'Y-dir')
    plt.plot(timeOsc,netForcePeriodBotUp1[0,:], label = 'X-dir')
    minValue = np.amin(netForcePeriodBotUp1)
    maxValue = np.amax(netForcePeriodBotUp1)
    plt.axis([0.0,timeFinal/10.0,-1.5*abs(minValue),1.5*abs(maxValue)])
    plt.legend(fontsize='x-small', loc='upper right')
    plt.savefig(dirName+'../BotUp1ForceRe'+str(Re)+'.png')
    plt.gcf().clear()
    plt.figure(1)
    plt.title('Net Force BotUp2: Re'+str(Re))
    plt.plot(timeOsc,netForcePeriodBotUp2[1,:],label = 'Y-dir')
    plt.plot(timeOsc,netForcePeriodBotUp2[0,:], label = 'X-dir')
    minValue = np.amin(netForcePeriodBotUp2)
    maxValue = np.amax(netForcePeriodBotUp2)
    plt.axis([0.0,timeFinal/10.0,-1.5*abs(minValue),1.5*abs(maxValue)])
    plt.legend(fontsize='x-small', loc='upper right')
    plt.savefig(dirName+'../BotUp2ForceRe'+str(Re)+'.png')
    plt.gcf().clear()
    plt.figure(2)
    plt.title('Net Force BotLow1: Re'+str(Re))
    plt.plot(timeOsc,netForcePeriodBotLow1[1,:],label = 'Y-dir')
    plt.plot(timeOsc,netForcePeriodBotLow1[0,:], label = 'X-dir')
    minValue = np.amin(netForcePeriodBotLow1)
    maxValue = np.amax(netForcePeriodBotLow1)
    plt.axis([0.0,timeFinal/10.0,-1.5*abs(minValue),1.5*abs(maxValue)])
    plt.legend(fontsize='x-small', loc='upper right')
    plt.savefig(dirName+'../BotLow1ForceRe'+str(Re)+'.png')
    plt.gcf().clear()
    plt.figure(3)
    plt.title('Net Force BotLow2: Re'+str(Re))
    plt.plot(timeOsc,netForcePeriodBotLow2[1,:],label = 'Y-dir')
    plt.plot(timeOsc,netForcePeriodBotLow2[0,:], label = 'X-dir')
    minValue = np.amin(netForcePeriodBotLow2)
    maxValue = np.amax(netForcePeriodBotLow2)
    plt.axis([0.0,timeFinal/10.0,-1.5*abs(minValue),1.5*abs(maxValue)])
    plt.legend(fontsize='x-small', loc='upper right')
    plt.savefig(dirName+'../BotLow2ForceRe'+str(Re)+'.png')
    plt.gcf().clear()
    #plt.show()
    
    #2) Find Force vector from last netForcePeriod value
    for i in range(2):                                           
        forceArr[0,0,i] = netForcePeriodBotUp1[i,-1]
        forceArr[0,1,i] = netForcePeriodBotUp2[i,-1]                                          
        forceArr[1,0,i] = netForcePeriodBotLow1[i,-1]
        forceArr[1,1,i] = netForcePeriodBotLow2[i,-1]
        print('BotUp1 Force[%i] = %f'%(i,forceArr[0,0,i]))
        print('BotUp2 Force[%i] = %f'%(i,forceArr[0,1,i]))
        print('BotLow1 Force[%i] = %f'%(i,forceArr[1,0,i]))
        print('BotLow2 Force[%i] = %f'%(i,forceArr[1,1,i]))

    return

def getFilePlots(nFPBU1,nFPBU2,nFPBL1,nFPBL2,Re):
    #Plot Net Force per Oscillation
    timeOsc = np.arange(0.0,timeFinal,PERIOD)
    plt.figure(0)
    plt.title('Net Force BotUp1: Re'+str(Re))
    plt.plot(timeOsc,nFPBU1[1,:],label = 'Y-dir')
    plt.plot(timeOsc,nFPBU1[0,:], label = 'X-dir')
    minValue = np.amin(nFPBU1)
    maxValue = np.amax(nFPBU1)
    plt.axis([0.0,timeFinal,-1.5*abs(minValue),1.5*abs(maxValue)])
    plt.legend(fontsize='x-small', loc='upper right')
    plt.savefig('BotUp1ForceRe'+str(Re)+'.png')
    plt.gcf().clear()
    plt.figure(1)
    plt.title('Net Force BotUp2: Re'+str(Re))
    plt.plot(timeOsc,nFPBU2[1,:],label = 'Y-dir')
    plt.plot(timeOsc,nFPBU2[0,:], label = 'X-dir')
    minValue = np.amin(nFPBU2)
    maxValue = np.amax(nFPBU2)
    plt.axis([0.0,timeFinal,-1.5*abs(minValue),1.5*abs(maxValue)])
    plt.legend(fontsize='x-small', loc='upper right')
    plt.savefig('BotUp2ForceRe'+str(Re)+'.png')
    plt.gcf().clear()
    plt.figure(2)
    plt.title('Net Force BotLow1: Re'+str(Re))
    plt.plot(timeOsc,nFPBL1[1,:],label = 'Y-dir')
    plt.plot(timeOsc,nFPBL1[0,:], label = 'X-dir')
    minValue = np.amin(nFPBL1)
    maxValue = np.amax(nFPBL1)
    plt.axis([0.0,timeFinal,-1.5*abs(minValue),1.5*abs(maxValue)])
    plt.legend(fontsize='x-small', loc='upper right')
    plt.savefig('BotLow1ForceRe'+str(Re)+'.png')
    plt.gcf().clear()
    plt.figure(3)
    plt.title('Net Force BotLow2: Re'+str(Re))
    plt.plot(timeOsc,nFPBL2[1,:],label = 'Y-dir')
    plt.plot(timeOsc,nFPBL2[0,:], label = 'X-dir')
    minValue = np.amin(nFPBL2)
    maxValue = np.amax(nFPBL2)
    plt.axis([0.0,timeFinal,-1.5*abs(minValue),1.5*abs(maxValue)])
    plt.legend(fontsize='x-small', loc='upper right')
    plt.savefig('BotLow2ForceRe'+str(Re)+'.png')
    plt.gcf().clear()
    #plt.show()
    return

def main3():
    #Make it more general! make it so that it can read in all the files from all the directories
    #1) Loop over all directories in list
    #2) Obtain hier_data files from each directory
    #3) Read files and calulcate net force
    #4) Create separate plot of net force for each case (like before) but now with COM!
    #5) End of loop: Compile all net forces for each directory in array
    #6) Plot force vectors for 1) Upper 2) Lower 3) COM
    
    nSims = 25
    nDimers = 2
    #Grid parameters
    xSize = 5
    ySize = 5   
    #Allocate Arrays
    #Directory List
    directoryList = [None]*nSims
    for i in range(xSize):
        for j in range(ySize):
            directoryList[5*i + j] = str(i)+str(j)
    print(directoryList)
    #Reynolds Number used
    Reynolds = [5,100]
    #Structure Names
    structureNames = ['botup','botlow']
    #Number of vertices in each structure
    nVert = np.zeros((xSize,ySize,len(structureNames),nDimers),dtype=int)
    #Net Force Over Period of All Bots
    #(xindx,yindx,structName,bot#,x or y)
    netForcePeriod = np.zeros((xSize,ySize,len(structureNames),nDimers,2)) #last period net force value in sim for each
    netForceCOMPeriod = np.zeros((xSize,ySize,2)) #Net force at COM for dimer2
    #COM of each structure
    comArr = np.zeros((xSize,ySize,len(structureNames),nDimers,2))
    
    #Loop over directories
    for Re in Reynolds:
        for dL in directoryList:
            dirName = '../'+dL+'/'
            print(dirName)
            #return
            xindx = int(dL[0])
            yindx = int(dL[1])
            #Obtain nvert for the structures used
            for i in range(len(structureNames)):
                for j in range(nDimers):
                    #Read vertex file line by line and save the values in a list using the \n delimiter
                    lines = [line.strip() for line in open(dirName+structureNames[i]+str(j+1)+'.vertex')]
                    #Break each list element into an array of numbers using the space delimiter
                    lines = [line.split() for line in lines]
                    nVert[xindx,yindx,i,j] = int(lines[0][0])
                    print('nvert[%i,%i,%i,%i] = %i'%(xindx,yindx,i,j,nVert[xindx,yindx,i,j]))
                    comArr[xindx,yindx,i,j,0] = float(lines[1][0])
                    comArr[xindx,yindx,i,j,1] = float(lines[1][1])
                    print('com[%i,%i,%i,%i,0] = %f'%(xindx,yindx,i,j,comArr[xindx,yindx,i,j,0]))
                    print('com[%i,%i,%i,%i,1] = %f'%(xindx,yindx,i,j,comArr[xindx,yindx,i,j,1]))                                 
            
            dirName = '../'+dL+'/hier_data_Re'+str(Re)+'/'
            getFilesForces(dirName,nVert[xindx,yindx,:,:],nDimers,Re,
                           comArr[xindx,yindx,:,:,:],netForcePeriod[xindx,yindx,:,:,:])
            
            #I have obtained the com and final net force for each dimer
            #Next I need to 1) plot the dimer pair in each directory separately with their net force vectors
            #2) Compile all of the forces and plot them
            #1)
            fig4 = plt.figure(4)
            ax4 = fig4.add_subplot(1,1,1,aspect=1)
            ax4.set_title('Steady State Force Direction: Re'+str(Re))
            #3) Part 1) Discs
            circle1 = plt.Circle((comArr[xindx,yindx,0,0,0],comArr[xindx,yindx,0,0,1]),0.15,color='r',fill=False)
            circle2 = plt.Circle((comArr[xindx,yindx,0,1,0],comArr[xindx,yindx,0,1,1]),0.15,color='b',fill=False)                                     
            circle3 = plt.Circle((comArr[xindx,yindx,1,0,0],comArr[xindx,yindx,1,0,1]),0.075,color='r',fill=False)
            circle4 = plt.Circle((comArr[xindx,yindx,1,1,0],comArr[xindx,yindx,1,1,1]),0.075,color='b',fill=False) 
            ax4.add_artist(circle1)
            ax4.add_artist(circle2)
            ax4.add_artist(circle3)
            ax4.add_artist(circle4)
            #3) Part 2) Force Arrows
            #Botup1
            ax4.arrow(comArr[xindx,yindx,0,0,0],comArr[xindx,yindx,0,0,1],
                      25.0*netForcePeriod[xindx,yindx,0,0,0],25.0*netForcePeriod[xindx,yindx,0,0,1],
                      head_width=0.025,head_length=0.05, fc='k',ec='k')
            #Botup2
            ax4.arrow(comArr[xindx,yindx,0,1,0],comArr[xindx,yindx,0,1,1],
                      25.0*netForcePeriod[xindx,yindx,0,1,0],25.0*netForcePeriod[xindx,yindx,0,1,1],
                      head_width=0.025,head_length=0.05, fc='k',ec='k')
            #BotLow1
            ax4.arrow(comArr[xindx,yindx,1,0,0],comArr[xindx,yindx,1,0,1],
                      25.0*netForcePeriod[xindx,yindx,1,0,0],25.0*netForcePeriod[xindx,yindx,1,0,1],
                      head_width=0.025,head_length=0.05, fc='k',ec='k')
            #BotLow2
            ax4.arrow(comArr[xindx,yindx,1,1,0],comArr[xindx,yindx,1,1,1],
                      25.0*netForcePeriod[xindx,yindx,1,1,0],25.0*netForcePeriod[xindx,yindx,1,1,1],
                      head_width=0.025,head_length=0.05, fc='k',ec='k')   
            ax4.axis([-0.5,1.2,-1.0,1.0])
            #ax4.axis('equal')
            fig4.savefig(dirName+'../SteadyStateRe'+str(Re)+'.png')
            fig4.clf()
            #plt.show()

        print('Completed reading all files! Time to make the compiled plots!')
        #return
        #Compile all of the forces and plot them
        fig5 = plt.figure(5)
        ax5 = fig5.add_subplot(1,1,1,aspect=1)
        ax5.set_title('Force Field Larger Disc: Re'+str(Re))
        fig6 = plt.figure(6)
        ax6 = fig6.add_subplot(1,1,1,aspect=1)
        ax6.set_title('Force Field Smaller Disc: Re'+str(Re))
        fig7 = plt.figure(7)
        ax7 = fig7.add_subplot(1,1,1,aspect=1)
        ax7.set_title('Force Field at COM: Re'+str(Re))
        #Plot Dimer 1
        BotUp1_5 = plt.Circle((comArr[0,0,0,0,0],comArr[0,0,0,0,1]),0.15,color='r',fill=False)
        BotLow1_5 = plt.Circle((comArr[0,0,1,0,0],comArr[0,0,1,0,1]),0.075,color='r',fill=False)
        BotUp1_6 = plt.Circle((comArr[0,0,0,0,0],comArr[0,0,0,0,1]),0.15,color='r',fill=False)
        BotLow1_6 = plt.Circle((comArr[0,0,1,0,0],comArr[0,0,1,0,1]),0.075,color='r',fill=False)
        BotUp1_7 = plt.Circle((comArr[0,0,0,0,0],comArr[0,0,0,0,1]),0.15,color='r',fill=False)
        BotLow1_7 = plt.Circle((comArr[0,0,1,0,0],comArr[0,0,1,0,1]),0.075,color='r',fill=False)
        ax5.add_artist(BotUp1_5)
        ax5.add_artist(BotLow1_5)
        ax6.add_artist(BotUp1_6)
        ax6.add_artist(BotLow1_6)
        ax7.add_artist(BotUp1_7)
        ax7.add_artist(BotLow1_7)
        for xindx in range(xSize):
            for yindx in range(ySize):
                #Plot dimer 2 forces
                #1) Larger disc forces
                ax5.arrow(comArr[xindx,yindx,0,1,0],comArr[xindx,yindx,0,1,1],
                          25.0*netForcePeriod[xindx,yindx,0,1,0],25.0*netForcePeriod[xindx,yindx,0,1,1],
                          head_width=0.025,head_length=0.05,fc='k',ec='k')
                #2) Smaller disc forces
                ax6.arrow(comArr[xindx,yindx,1,1,0],comArr[xindx,yindx,1,1,1],
                          25.0*netForcePeriod[xindx,yindx,1,1,0],25.0*netForcePeriod[xindx,yindx,1,1,1],
                          head_width=0.025,head_length=0.05,fc='k',ec='k')
                #3) COM forces
                netForceCOMPeriod[xindx,yindx,0] = netForcePeriod[xindx,yindx,0,1,0] + netForcePeriod[xindx,yindx,1,1,0]
                netForceCOMPeriod[xindx,yindx,1] = netForcePeriod[xindx,yindx,0,1,1] + netForcePeriod[xindx,yindx,1,1,1]
                ax7.arrow(comArr[xindx,yindx,0,1,0],comArr[xindx,yindx,0,1,1],
                          25.0*netForceCOMPeriod[xindx,yindx,0],25.0*netForceCOMPeriod[xindx,yindx,1],
                          head_width=0.025,head_length=0.05,fc='k',ec='k')
                print('(%i,%i)' %(xindx,yindx))
    
        ax5.axis([-0.5,1.2,-0.7,1.0])
        ax6.axis([-0.5,1.2,-0.7,1.0])
        ax7.axis([-0.5,1.2,-0.7,1.0])
        fig5.savefig('ForceFieldLargeRe'+str(Re)+'.png')
        fig6.savefig('ForceFieldSmallRe'+str(Re)+'.png')
        fig7.savefig('ForceFieldCOMRe'+str(Re)+'.png')
        fig5.clf()
        fig6.clf()
        fig7.clf()
        #return
        
    return

def reject_outliers(data):
    m = 2
    return data[abs(data - np.mean(data)) < m * np.std(data)]

def get_distance_bw(Arr1,Arr2):
    distSquared = (Arr1[0,:]-Arr2[0,:])*(Arr1[0,:]-Arr2[0,:]) + (Arr1[1,:] - Arr2[1,:])*(Arr1[1,:]-Arr2[1,:])
    return np.sqrt(distSquared)

def main():    
    
    #Data Dump Number
    fileNumber = np.arange(100,500100,100)
    print("fileNumber[10] = ",fileNumber[10])
    
    #Constants
    dt = 1.0e-5
    Ks = "Ks1e6"
    Re = "Re5"
    
    #Allocate Arrays for pos, time, and vel
    posArr = np.zeros((2,2,len(fileNumber)+1)) #struct,axis,line
    topArr = np.zeros((2,2,len(fileNumber)+1))
    botArr = np.zeros((2,2,len(fileNumber)+1))
    timeArr = np.zeros(len(fileNumber)+1)
    velArr = np.zeros((2,2,len(fileNumber))) #struct,axis,line
    magVelArr = np.zeros((2,len(fileNumber))) #struct,line
    desiredLength = np.zeros(len(fileNumber)+1)
    currentLength = np.zeros(len(fileNumber)+1)
    slopeArr = np.zeros((3,len(fileNumber)+1)) #0 = upper top to upper center 1 = lower top to lower center 2 = upper center to lower center
    
    #Initial Conditions
    posArr[0,0,0] = 0.0
    posArr[0,1,0] = 0.0
    posArr[1,0,0] = 0.0
    posArr[1,1,0] = -0.375
    topArr[0,0,0] = 0.0
    topArr[0,1,0] = 0.15
    topArr[1,0,0] = 0.0
    topArr[1,1,0] = -0.30
    botArr[0,0,0] = 0.0
    botArr[0,1,0] = -0.15
    botArr[1,0,0] = 0.0
    botArr[1,1,0] = -0.45
    timeArr[0] = 0.0
    desiredLength[0] = 0.375
    currentLength[0] = 0.375
    
    #Determine Steady Velocity of disks
    #steadyVel = np.zeros(len(Reynolds))
    fileCount = 0
    
    plot1Title = "Disk deltaY vs Time: "+str(Re)+" "+str(Ks)
    plot2Title = "Desired/Current Length vs Time: "+str(Re)+" "+str(Ks)
    plot3Title = "Percent Difference Spring Length vs. Time: "+str(Re)+" "+str(Ks)
    
    for fN in fileNumber:
        
        #Iterate fileCount
        fileCount += 1
        print("fN = ",fN)
    
        #pwd including file name
        pwd = "/Users/thomas/sfw/visitFiles/thomas/1bot/StandardIB/RefinementStudy/N256/0.5dX"
        if(fN < 1000):
            PDpwd = pwd+"/positionData"+str(Ks)+"t5"+str(Re)+"/X.00"+str(fN)
        elif(fN >= 1000 and fN < 10000):
            PDpwd = pwd+"/positionData"+str(Ks)+"t5"+str(Re)+"/X.0"+str(fN)
        else:
            PDpwd = pwd+"/positionData"+str(Ks)+"t5"+str(Re)+"/X."+str(fN)
    
        #Read pd.txt file line by line and save the values in a list using the \n delimiter
        lines = [line.strip() for line in open(PDpwd)]
        #Break each list element into an array of numbers using the space delimiter
        lines = [line.split() for line in lines]
        
        #Store COM positions and time for upper and lower discs
        #Upper disc
        posArr[0,0,fileCount] = float(lines[7][0])
        posArr[0,1,fileCount] = float(lines[8][0])
        #print("ytop = ",posArr[0,1,fileCount])
        topArr[0,0,fileCount] = float(lines[3][0])
        topArr[0,1,fileCount] = float(lines[4][0])
        botArr[0,0,fileCount] = float(lines[11][0])
        botArr[0,1,fileCount] = float(lines[12][0])
        #Lower disc
        posArr[1,0,fileCount] = float(lines[17][0])
        posArr[1,1,fileCount] = float(lines[18][0])
        topArr[1,0,fileCount] = float(lines[13][0])
        topArr[1,1,fileCount] = float(lines[14][0])
        botArr[1,0,fileCount] = float(lines[21][0])
        botArr[1,1,fileCount] = float(lines[22][0])
        #Time
        timeArr[fileCount] = (dt*fN) 

    #End of file extraction
        
    for i in range(1,len(fileNumber)+1):
        #Current Spring Length
        xDist = posArr[0,0,i] - posArr[1,0,i]
        yDist = posArr[0,1,i] - posArr[1,1,i]
        currentLength[i] = np.sqrt(xDist*xDist + yDist*yDist)
        #Desired Spring Length
        desiredLength[i] = 0.375 + 0.5*0.3*0.8*np.sin(2.0*np.pi*8.0*timeArr[i])
        #Calculate slopes formed by discs' axes and displacement b/w centers
        slopeArr[0,i] = (topArr[0,0,i] - posArr[0,0,i])/(topArr[0,1,i] - posArr[0,1,i])
        slopeArr[1,i] = (topArr[1,0,i] - posArr[1,0,i])/(topArr[1,1,i] - posArr[1,1,i])
        slopeArr[2,i] = (posArr[0,0,i] - posArr[1,0,i])/(posArr[0,1,i] - posArr[1,1,i])   
        if(i > 0 and i < len(fileNumber)):
            #Vel of Struct 0
            velArr[0,0,i] = (posArr[0,0,i] - posArr[0,0,i-1])/(timeArr[i] - timeArr[i-1])
            velArr[0,1,i] = (posArr[0,1,i] - posArr[0,1,i-1])/(timeArr[i] - timeArr[i-1])
            magVelArr[0,i] = np.sqrt(velArr[0,0,i]*velArr[0,0,i] + velArr[0,1,i]*velArr[0,1,i])
            #Vel of Struct 1
            velArr[1,0,i] = (posArr[1,0,i] - posArr[1,0,i-1])/(timeArr[i] - timeArr[i-1])
            velArr[1,1,i] = (posArr[1,1,i] - posArr[1,1,i-1])/(timeArr[i] - timeArr[i-1])
            magVelArr[1,i] = np.sqrt(velArr[1,0,i]*velArr[1,0,i] + velArr[1,1,i]*velArr[1,1,i])
    
        #Determine a final velocity in the simulation using linear regression
        #Find slope of trajectory for last second of simulation
        '''startTimeIndex = len(lines) - 10001
        storeMaxPos = []
        storeMaxTime = []
        index = []
        #Attempt 2: Sample out maximums, form line between them
        for ii in range(startTimeIndex,len(lines)):
            #Check Velocity
            if(np.sign(velArr[1,1,ii]) != np.sign(velArr[1,1,ii-1]) and velArr[1,1,ii] < 0.0):
                #There is a maximum at velArr[1,1,ii-1]
                storeMaxTime.append(timeArr[ii-1])
                storeMaxPos.append(posArr[1,1,ii-1])
        #Filter out Outliers
        m = 2
        mean = np.mean(storeMaxPos)
        std = np.std(storeMaxPos)
        filteredPos = [e for e in storeMaxPos if (mean - m*std < e < mean + m*std)]
        outlierPos = [e for e in storeMaxPos if (mean - m*std > e or e > mean + m*std)]
        for i in range(len(outlierPos)):
            for j in range(len(storeMaxPos)):
                if(storeMaxPos[j] == outlierPos[i]):
                    #Store Index of outlier
                    index.append(j)
                
        print("filtered = ",filteredPos)
        print("outlier = ",outlierPos)
        print("raw = ",storeMaxPos)
        for j in range(len(index)):
            #Remove outlier from stored Pos and Time values
            storeMaxTime.remove(storeMaxTime[index[j]])
            storeMaxPos.remove(storeMaxPos[index[j]])
        print("new raw = ",storeMaxPos)
        print("new time = ",storeMaxTime)
        #Plot Max Pos vs time
        plt.figure(6)
        plt.plot(storeMaxTime,storeMaxPos,'ro')
        plt.show()
        plt.gcf().clear()
        #Find Linear Regression of the Max Pos vs time
        slope, intercept, r_value, p_value, std_err = stats.linregress(storeMaxTime,storeMaxPos)
        print("r_value for Re"+str(Re)+" = ",r_value)
        steadyVel[ReCount] = slope'''
    
    #Plot ypos vs time
    plt.figure(0)
    plt.title("Disk Y-Pos vs. Time: "+str(Re)+" "+str(Ks))
    #disk1Line = slope*timeArr[startTimeIndex:] + intercept
    #plt.plot(timeArr[startTimeIndex:],disk1Line, 'r', label = "Steady Vel")
    plt.plot(timeArr,[-0.375]*(len(fileNumber)+1),'k--')
    plt.plot(timeArr,posArr[0,1,:])
    plt.plot(timeArr,[0.0]*(len(fileNumber)+1),'k--')
    plt.plot(timeArr,posArr[1,1,:])
    plt.xlabel("time (s)")
    plt.ylabel("Y-pos (m)")
    #plt.axis([0.0,timeArr[-1],np.minimum(posArr[1,1,-1] - 0.5,-1.0),np.maximum(posArr[0,1,-1] + 0.5,1.0)])
    plt.axis([0.0,5.0,-1.0,0.1])
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.savefig(pwd+"/diskTrajectory"+str(Re)+str(Ks)+".png")
    plt.gcf().clear()
    '''#Plot deltaY vs time
    plt.figure(1)
    plt.title(plot1Title)
    plt.plot(timeArr,posArr[0,1,:],label = "disk0")
    plt.plot(timeArr,posArr[1,1,:]+0.375,label = "disk1")
    plt.plot(timeArr,posArr[0,1,:] - posArr[1,1,:] - 0.375, label = "surf-surf")
    plt.xlabel("time (s)")
    plt.ylabel("deltaY (m)")
    plt.axis([4.0,5.0,-0.15,0.15])
    plt.legend(loc="lower left")
    plt.tight_layout()
    plt.savefig(pwd+"/deltaY"+str(Re)+str(Ks)+".png")
    plt.gcf().clear()
    #Plot Desired/Current vs time
    plt.figure(2)
    plt.title(plot2Title)
    plt.plot(timeArr,desiredLength-0.375,'r', label = "desired")
    plt.plot(timeArr,currentLength-0.375,'b', label = "current")
    plt.xlabel("time (s)")
    plt.ylabel("Spring Length (m)")
    plt.axis([4.0,5.0,-0.15,0.15])
    plt.legend(loc="lower left")
    plt.tight_layout()
    plt.savefig(pwd+"/springLength"+str(Re)+str(Ks)+".png")
    plt.gcf().clear()'''
    #Plot Percent Diff for oscillation (desired vs current)
    plt.figure(3)
    plt.title(plot3Title)
    percentDiff = 1.0E2*(desiredLength - currentLength)/desiredLength
    plt.plot(timeArr,percentDiff)
    plt.xlabel("time (s)")
    plt.ylabel("%Diff")
    plt.axis([0.0,1.0,-1.0,0.25])
    plt.tight_layout()
    plt.savefig(pwd+"/percentDiff"+str(Re)+str(Ks)+".png")
    plt.gcf().clear()
    #Plot X-dist vs time
    plt.figure(4)
    plt.title("X-Dist vs time: "+str(Re)+" "+str(Ks))
    plt.plot(timeArr,posArr[0,0,:], label = "upper")
    plt.plot(timeArr,posArr[1,0,:], label = "lower")
    plt.xlabel("time (s)")
    plt.ylabel("X-Dist (m)")
    plt.axis([0.0,5.0,-0.05,0.05])
    plt.tight_layout()
    plt.savefig(pwd+"/X-DistvTime"+str(Re)+str(Ks)+".png")
    plt.gcf().clear()
    '''#Plot slope vs time
    plt.figure(5)
    plt.title("Slopes UT->UC vs time: "+str(Re)+" "+str(Ks))
    plt.plot(timeArr,slopeArr[0,:], label = "ut->uc")
    plt.plot(timeArr,slopeArr[1,:], label = "lt->lc")
    plt.plot(timeArr,slopeArr[2,:], label = "uc->lc")
    plt.xlabel("time (s)")
    plt.ylabel("Slope ")
    plt.axis([4.0,5.0,-0.05,0.05])
    plt.legend(fontsize="x-small")
    plt.tight_layout()
    plt.savefig(pwd+"/SlopevTime"+str(Re)+str(Ks)+".png")
    plt.gcf().clear()
    #Plot Slope Percent Differences
    plt.figure(6)
    plt.subplot(211)
    plt.title("Angle Diff UTUC->UCLC vs time: "+str(Re)+" "+str(Ks))
    #slopePercentDiff = 1.0E2*(slopeArr[2,:] - slopeArr[0,:])/slopeArr[2,:]
    slopeAngle = (topArr[0,0,:] - posArr[0,0,:])*(posArr[0,0,:] - posArr[1,0,:]) + (topArr[0,1,:] - posArr[0,1,:])*(posArr[0,1,:] - posArr[1,1,:])
    DistCC = get_distance_bw(posArr[0,:,:],posArr[1,:,:])
    DistTC = get_distance_bw(topArr[0,:,:],posArr[0,:,:])
    CosAngle = np.zeros(len(timeArr))
    InvCos = np.zeros(len(timeArr))
    for i in range(len(timeArr)):
        CosAngle[i] = slopeAngle[i]/(DistCC[i]*DistTC[i])
        print("CosAngle = ",CosAngle[i])
        InvCos[i] = np.arccos(CosAngle[i])
        print("InvCos = ",InvCos[i])
    plt.plot(timeArr,InvCos*180.0/np.pi)
    plt.xlabel("time (s)")
    plt.ylabel("theta")
    plt.axis([0.0,1.0,-0.01,0.01])
    plt.subplot(212)
    plt.title("Angle Diff LTLC->UCLC vs time: "+str(Re)+" "+str(Ks))
    #slopePercentDiff = 1.0E2*(slopeArr[2,:] - slopeArr[1,:])/slopeArr[2,:]
    slopeAngle = (topArr[1,0,:] - posArr[1,0,:])*(posArr[0,0,:] - posArr[1,0,:]) + (topArr[1,1,:] - posArr[1,1,:])*(posArr[0,1,:] - posArr[1,1,:])
    print("dotProduct = ",slopeAngle[0])
    print("dotProduct = ",slopeAngle[100])
    DistCC = get_distance_bw(posArr[0,:,:],posArr[1,:,:])
    print("DistCC = ",DistCC[0])
    print("DistCC = ",DistCC[100])
    DistTC = get_distance_bw(topArr[1,:,:],posArr[1,:,:])
    print("DistTC = ",DistTC[0])
    print("DistTC = ",DistTC[100])
    for i in range(len(timeArr)):
        CosAngle[i] = slopeAngle[i]/(DistCC[i]*DistTC[i])
        InvCos[i] = np.arccos(CosAngle[i])
        #print("cos(theta) = ",CosAngle[i])
    plt.plot(timeArr,InvCos*180.0/np.pi)
    plt.xlabel("time (s)")
    plt.ylabel("theta ")
    plt.axis([0.0,1.0,-0.01,0.01])
    plt.tight_layout()
    plt.savefig(pwd+"/SlopePercentDiff"+str(Re)+str(Ks)+".png")
    plt.gcf().clear()'''
    #X-Dist Change (Top)
    plt.figure(7)
    plt.title("X-Dist UT/LT vs time: "+str(Re)+" "+str(Ks))
    Xdist = topArr[0,0,:] - posArr[0,0,:]
    plt.plot(timeArr,1.0e2*Xdist/0.15, label = "Upper")
    Xdist = topArr[1,0,:] - posArr[1,0,:]
    plt.plot(timeArr,1.0e2*Xdist/0.075, label = "lower")
    plt.xlabel("time (s)")
    plt.ylabel("X-Dist (m)")
    plt.axis([4.0,5.0,-10.0,10.0])
    plt.legend(fontsize='x-small')
    plt.tight_layout()
    plt.savefig(pwd+"/X-DistTopvTime"+str(Re)+str(Ks)+".png")
    plt.gcf().clear()
    #Y-Dist Change (Top)
    plt.figure(8)
    plt.title("Y-Dist UT/LT vs time: "+str(Re)+" "+str(Ks))
    Ydist = topArr[0,1,:] - posArr[0,1,:] - 0.15
    plt.plot(timeArr,1.0e2*Ydist/0.15, label = "Upper")
    Ydist = topArr[1,1,:] - posArr[1,1,:] - 0.075
    plt.plot(timeArr,1.0e2*Ydist/0.075, label = "lower")
    plt.xlabel("time (s)")
    plt.ylabel("Y-Dist (m)")
    plt.axis([4.0,5.0,-0.05,0.05])
    plt.legend(fontsize='x-small')
    plt.tight_layout()
    plt.savefig(pwd+"/Y-DistTopvTime"+str(Re)+str(Ks)+".png")
    plt.gcf().clear()
    plt.figure(9)
    #Rotation on top and bottom discs
    plt.title("CounterClockwise (+) or Clockwise(-) vs. Time: "+str(Re)+" "+str(Ks))
    SinAngle = np.zeros(len(timeArr))
    InvSin = np.zeros(len(timeArr))
    DistCC = get_distance_bw(posArr[0,:,:],posArr[1,:,:])
    for i in range(len(timeArr)):
        SinAngle[i] = (posArr[0,0,i] - posArr[1,0,i])/DistCC[i]
        print("SinAngle = ",SinAngle[i])
        InvSin[i] = np.arccos(SinAngle[i])
        print("InvSin = ",InvSin[i])
    dInvSin = np.zeros(len(timeArr))
    dInvSin[0] = 0.0
    for i in range(1,len(timeArr)):
        dInvSin[i] = (InvSin[i] - InvSin[i-1])/(1.0e-3)
    plt.plot(timeArr,dInvSin,label = "rotation")
    plt.plot(timeArr,(posArr[1,0,:]-posArr[0,0,:])*100.0,label="XDist")
    plt.xlabel("time (s)")
    plt.ylabel("Rate of Rotation (rad/s)")
    plt.axis([0.0,1.0,-0.5,0.5])
    plt.tight_layout()
    plt.legend(fontsize="x-small")
    plt.savefig(pwd+"/RotationvTime"+str(Re)+str(Ks)+".png")
    plt.gcf().clear()
      
#-------------_END_PROGRAM_-----------------
#main()
#main2()        
main3()
