#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 14:31:19 2018

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
    dumpInterval = 500
    maxFileNumber = int(timeFinal/(dumpInterval*dt))*dumpInterval
    fileNumber = np.arange(500,maxFileNumber+dumpInterval,dumpInterval)
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
        netForcePeriodBotUp1[0,i] = np.sum(netForceBotUp1[0,25*i:25*(i+1)])/(PERIOD/(dumpInterval*dt))
        netForcePeriodBotUp2[0,i] = np.sum(netForceBotUp2[0,25*i:25*(i+1)])/(PERIOD/(dumpInterval*dt))
        netForcePeriodBotLow1[0,i] = np.sum(netForceBotLow1[0,25*i:25*(i+1)])/(PERIOD/(dumpInterval*dt))
        netForcePeriodBotLow2[0,i] = np.sum(netForceBotLow2[0,25*i:25*(i+1)])/(PERIOD/(dumpInterval*dt))
        netForcePeriodBotUp1[1,i] = np.sum(netForceBotUp1[1,25*i:25*(i+1)])/(PERIOD/(dumpInterval*dt))
        netForcePeriodBotUp2[1,i] = np.sum(netForceBotUp2[1,25*i:25*(i+1)])/(PERIOD/(dumpInterval*dt))
        netForcePeriodBotLow1[1,i] = np.sum(netForceBotLow1[1,25*i:25*(i+1)])/(PERIOD/(dumpInterval*dt))
        netForcePeriodBotLow2[1,i] = np.sum(netForceBotLow2[1,25*i:25*(i+1)])/(PERIOD/(dumpInterval*dt))

    #Plot Net Force per Oscillation
    timeOsc = np.arange(0.0,timeFinal,PERIOD)
    plt.figure(0)
    plt.title('Net Force BotUp1: Re'+str(Re))
    plt.plot(timeOsc,netForcePeriodBotUp1[1,:],label = 'Y-dir')
    plt.plot(timeOsc,netForcePeriodBotUp1[0,:], label = 'X-dir')
    minValue = np.amin(netForcePeriodBotUp1)
    maxValue = np.amax(netForcePeriodBotUp1)
    plt.axis([0.0,timeFinal,-1.5*abs(minValue),1.5*abs(maxValue)])
    plt.legend(fontsize='x-small', loc='upper right')
    plt.savefig(dirName+'../BotUp1ForceRe'+str(Re)+'.png')
    plt.gcf().clear()
    plt.figure(1)
    plt.title('Net Force BotUp2: Re'+str(Re))
    plt.plot(timeOsc,netForcePeriodBotUp2[1,:],label = 'Y-dir')
    plt.plot(timeOsc,netForcePeriodBotUp2[0,:], label = 'X-dir')
    minValue = np.amin(netForcePeriodBotUp2)
    maxValue = np.amax(netForcePeriodBotUp2)
    plt.axis([0.0,timeFinal,-1.5*abs(minValue),1.5*abs(maxValue)])
    plt.legend(fontsize='x-small', loc='upper right')
    plt.savefig(dirName+'../BotUp2ForceRe'+str(Re)+'.png')
    plt.gcf().clear()
    plt.figure(2)
    plt.title('Net Force BotLow1: Re'+str(Re))
    plt.plot(timeOsc,netForcePeriodBotLow1[1,:],label = 'Y-dir')
    plt.plot(timeOsc,netForcePeriodBotLow1[0,:], label = 'X-dir')
    minValue = np.amin(netForcePeriodBotLow1)
    maxValue = np.amax(netForcePeriodBotLow1)
    plt.axis([0.0,timeFinal,-1.5*abs(minValue),1.5*abs(maxValue)])
    plt.legend(fontsize='x-small', loc='upper right')
    plt.savefig(dirName+'../BotLow1ForceRe'+str(Re)+'.png')
    plt.gcf().clear()
    plt.figure(3)
    plt.title('Net Force BotLow2: Re'+str(Re))
    plt.plot(timeOsc,netForcePeriodBotLow2[1,:],label = 'Y-dir')
    plt.plot(timeOsc,netForcePeriodBotLow2[0,:], label = 'X-dir')
    minValue = np.amin(netForcePeriodBotLow2)
    maxValue = np.amax(netForcePeriodBotLow2)
    plt.axis([0.0,timeFinal,-1.5*abs(minValue),1.5*abs(maxValue)])
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

def main():
    #Make it more general! make it so that it can read in all the files from all the directories
    #1) Loop over all directories in list
    #2) Obtain hier_data files from each directory
    #3) Read files and calulcate net force
    #4) Create separate plot of net force for each case (like before) but now with COM!
    #5) End of loop: Compile all net forces for each directory in array
    #6) Plot force vectors for 1) Upper 2) Lower 3) COM
    
    nDimers = 2
    #Grid parameters
    xSize = 3
    ySize = 3  
    nSims = int(xSize*ySize)
    #Allocate Arrays
    #Dimer Location
    dimerLocation = ['Above','Below']
    #Directory List
    directoryList = [None]*nSims
    for i in range(xSize):
        for j in range(ySize):
            directoryList[xSize*i + j] = str(i)+str(j)
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
    
    #File to save all Data
    f = open('../'+dimerLocation[0]+'/'+dimerLocation[0]+'DataFile.txt','w')
    g = open('../'+dimerLocation[1]+'/'+dimerLocation[1]+'DataFile.txt','w')
    #Format:        xpos ypos Fx Fy
    #      Re5
    #      botup1
    #      botlow1
    #      botup2
    #      botlow2
    #      bot2COM
    #      Repeat
    #      Re100 Repeat
    
    #Loop over directories
    for Re in Reynolds:
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
        for dLoc in dimerLocation:
            for dL in directoryList:
                dirName = '../'+dLoc+'/'+dL+'/'+'Re'+str(Re)+'/'
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
                
                dirName = '../'+dLoc+'/'+dL+'/Re'+str(Re)+'/hier_data_Re'+str(Re)+'/'
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
                ax4.axis([0.0,1.0,-1.0,1.5])
                #ax4.axis('equal')
                fig4.savefig(dirName+'../SteadyStateRe'+str(Re)+'.png')
                fig4.clf()
                #plt.show()

            print('Completed reading all files! Time to make the compiled plots!')
            #return
        
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

                    #Write to Data File
                    if(dLoc == 'Above'):
                        f.write('%.3f %.3f %.3f %.3f\n'%(comArr[xindx,yindx,0,0,0],comArr[xindx,yindx,0,0,1],
                                                       netForcePeriod[xindx,yindx,0,0,0],netForcePeriod[xindx,yindx,0,0,1]))
                        f.write('%.3f %.3f %.3f %.3f\n'%(comArr[xindx,yindx,1,0,0],comArr[xindx,yindx,1,0,1],
                                                       netForcePeriod[xindx,yindx,1,0,0],netForcePeriod[xindx,yindx,1,0,1]))
                        f.write('%.3f %.3f %.3f %.3f\n'%(comArr[xindx,yindx,0,1,0],comArr[xindx,yindx,0,1,1],
                                                       netForcePeriod[xindx,yindx,0,1,0],netForcePeriod[xindx,yindx,0,1,1]))
                        f.write('%.3f %.3f %.3f %.3f\n'%(comArr[xindx,yindx,1,1,0],comArr[xindx,yindx,1,1,1],
                                                       netForcePeriod[xindx,yindx,1,0,0],netForcePeriod[xindx,yindx,1,0,1]))
                        f.write('%.3f %.3f %.3f %.3f\n'%(comArr[xindx,yindx,0,1,0],comArr[xindx,yindx,0,1,1]-0.063,
                                                       netForceCOMPeriod[xindx,yindx,0],netForceCOMPeriod[xindx,yindx,1]))
                    if(dLoc == 'Below'):
                        g.write('%.3f %.3f %.3f %.3f\n'%(comArr[xindx,yindx,0,0,0],comArr[xindx,yindx,0,0,1],
                                                       netForcePeriod[xindx,yindx,0,0,0],netForcePeriod[xindx,yindx,0,0,1]))
                        g.write('%.3f %.3f %.3f %.3f\n'%(comArr[xindx,yindx,1,0,0],comArr[xindx,yindx,1,0,1],
                                                       netForcePeriod[xindx,yindx,1,0,0],netForcePeriod[xindx,yindx,1,0,1]))
                        g.write('%.3f %.3f %.3f %.3f\n'%(comArr[xindx,yindx,0,1,0],comArr[xindx,yindx,0,1,1],
                                                       netForcePeriod[xindx,yindx,0,1,0],netForcePeriod[xindx,yindx,0,1,1]))
                        g.write('%.3f %.3f %.3f %.3f\n'%(comArr[xindx,yindx,1,1,0],comArr[xindx,yindx,1,1,1],
                                                       netForcePeriod[xindx,yindx,1,0,0],netForcePeriod[xindx,yindx,1,0,1]))
                        g.write('%.3f %.3f %.3f %.3f\n'%(comArr[xindx,yindx,0,0,0],comArr[xindx,yindx,0,0,1]-0.063,
                                                       netForceCOMPeriod[xindx,yindx,0],netForceCOMPeriod[xindx,yindx,1]))
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
    
        ax5.axis([0.0,1.0,-1.0,1.5])
        ax6.axis([0.0,1.0,-1.0,1.5])
        ax7.axis([0.0,1.0,-1.0,1.5])
        fig5.savefig('../ForceFieldLargeABRe'+str(Re)+'.png')
        fig6.savefig('../ForceFieldSmallABRe'+str(Re)+'.png')
        fig7.savefig('../ForceFieldCOMABRe'+str(Re)+'.png')
        fig5.clf()
        fig6.clf()
        fig7.clf()        
        
        #return
    f.close()
    g.close()
        
    return
      
#-------------_END_PROGRAM_-----------------        
main()
