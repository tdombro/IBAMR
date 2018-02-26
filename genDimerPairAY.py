#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 12:05:02 2018

@author: thomas
"""

#1)Creates an array of spatial distances
#   2 dimers (Fixed positions, but distance b/w them changes)
#2) Random Orientation (Not yet, keep vertical atm) 
#3) zips all files (.vertex,.spring,.beam,.target)

import numpy as np
import matplotlib.pyplot as plt
import os
import zipfile
from shutil import copyfile

def DimerLocation(Nbots,structureName,COM,dirName,Re,ax0):
    #First dimer always located at COM = [0.5,0.5]
    #Sets up a grid of locations for the second dimer
    #Does not include positions directly above and below the first dimer
    
    #Loops to create grid and make directories
    for i in range(len(structureName)):    
        #Determine # vertices in current structure
        #Read vertex file line by line and save the values in a list using the \n delimiter
        lines = [line.strip() for line in open(structureName[i]+str(1)+'.vertex')]
        #Break each list element into an array of numbers using the space delimiter
        lines = [line.split() for line in lines]
        nvert = int(lines[0][0])
        print('nvert = ',nvert)
        
        #Allocate Arrays
        arrPosition = np.zeros((Nbots,2,nvert))
        
        #Store vertex coord values for first dimer
        vertexNumber = 1
        for k in range(1,nvert+1):
            arrPosition[vertexNumber-1,0,k-1] = float(lines[k][0])
            arrPosition[vertexNumber-1,1,k-1] = float(lines[k][1])
            #Translate Vertex coordinates to fit in current grid cell
            arrPosition[vertexNumber-1,1,k-1] -= COM[0,1]
            if(Re == '5'):
                arrPosition[vertexNumber-1,0,k-1] += -0.305
            if(Re == '100'):
                arrPosition[vertexNumber-1,0,k-1] += COM[0,1]
            
        #Store vertex coord values for second dimer
        vertexNumber += 1
        for k in range(1,nvert+1):
            arrPosition[vertexNumber-1,0,k-1] = float(lines[k][0])
            arrPosition[vertexNumber-1,1,k-1] = float(lines[k][1])
            #Translate Vertex coordinates to fit in current grid cell
            if(Re == '5'):
                arrPosition[vertexNumber-1,0,k-1] += 2.0 - np.sqrt(3.0) - 0.1525*np.sqrt(3.0)
                arrPosition[vertexNumber-1,1,k-1] += 1.0 + 0.375 - 0.1525
            if(Re == '100'):
                arrPosition[vertexNumber-1,0,k-1] += 2.0 - 1.45*np.sqrt(3.0) + COM[0,1]*np.sqrt(3.0)/2.0
                arrPosition[vertexNumber-1,1,k-1] += 1.45 - 3.0*COM[0,1]/2.0
            
        #Create Zip file of all .vertex, .spring, and .target files of interest
        #Copy .spring and .target into directory of interest
        for j in range(1,Nbots+1):
            copyfile(structureName[i]+str(1)+'.spring',dirName+structureName[i]+str(j)+'.spring')
            if(structureName[i] == 'skeleton'):
                copyfile(structureName[i]+str(1)+'.beam',dirName+structureName[i]+str(j)+'.beam')

        #Plot Skeleton
        ax0.plot(arrPosition[0,0,:],arrPosition[0,1,:],'ro')
        ax0.plot(arrPosition[1,0,:],arrPosition[1,1,:],'ro')
        ax0.axis([-1.0,1.0,0.0,2.0])
        ax0.axis('equal')
    
    #Translate COM (only for dimer 2)
    print('b4: COMxbot1 = ',COM[0,0])
    print('b4: COMx = ',COM[1,0])
    print('b4: COMy = ',COM[1,1])
    if(Re == '5'):
        COM[0,0] += -0.305
        COM[0,1] = 0.0
        COM[1,0] += 2.0 - np.sqrt(3.0) - 0.1525*np.sqrt(3.0)
        COM[1,1] += 1.0 + 0.375 - 0.1525
    if(Re == '100'):
        COM[1,0] += 2.0 - 1.45*np.sqrt(3.0) + COM[0,1]*np.sqrt(3.0)/2.0
        COM[1,1] += 1.45 - 3.0*COM[0,1]/2.0
        COM[0,0] += COM[0,1]
        COM[0,1] = 0.0
    print('a4: COMx = ',COM[1,0])
    print('a4: COMy = ',COM[1,1])
    
    return (arrPosition,nvert,COM)           

def RotateVertices(Nbots,nvert,structureName,dirName,gridPosition,COM,Re):  
    #Rotates the set up grid by a randomly generated angle about the COM
    
    #Plot Rotated Dimers
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(1,1,1,aspect=1)
    ax1.set_title('Rotated Dimers Re'+str(Re))
    
    #Allocate Arrays
    rotationMatrix = np.zeros((Nbots,2,2))
    theta = np.zeros(Nbots)
    rotatedPosition = np.zeros((Nbots,2,nvert))
    #COM = np.zeros((Nbots,2))
    
    #Generate Angles
    if(Re == '5'):
        theta[0] = np.pi/2.0
        theta[1] = np.pi/3.0
    if(Re == '100'):
        theta[0] = np.pi/2.0 - np.pi
        theta[1] = np.pi/3.0 - np.pi
    for i in range(Nbots):
        rotationMatrix[i,0,0] = np.cos(theta[i])
        rotationMatrix[i,0,1] = -1.0*np.sin(theta[i])
        rotationMatrix[i,1,0] = np.sin(theta[i])
        rotationMatrix[i,1,1] = np.cos(theta[i])
    
    '''#Generate Random Angles
    np.random.seed(seed)
    for i in range(Nbots):
        theta[i] = 2.0*np.pi*np.random.rand(1)
        print('theta[%i] = %.3e' %(i,theta[i]))
        rotationMatrix[i,0,0] = np.cos(theta[i])
        rotationMatrix[i,0,1] = -1.0*np.sin(theta[i])
        rotationMatrix[i,1,0] = np.sin(theta[i])
        rotationMatrix[i,1,1] = np.cos(theta[i])'''     

    #Rotate Grid Structures around COM based off of above angle 
    for i in range(int(Nbots)):
        #open newly rotated vertex file
        f = open(dirName+structureName+str(i+1)+'.vertex','w')
        f.write('%i\n'%nvert)
        
        #Rotate coord values about center of mass
        for j in range(nvert):
            rotatedPosition[i,:,j] = rotationMatrix[i].dot(gridPosition[i,:,j] - COM[i,:])
            rotatedPosition[i,:,j] += COM[i,:]
            if(j == nvert - 1):
                f.write('%.5e %.5e'%(rotatedPosition[i,0,j],rotatedPosition[i,1,j]))
            else:
                f.write('%.5e %.5e\n' %(rotatedPosition[i,0,j],rotatedPosition[i,1,j]))
        f.close()
        
        #Plot Rotated Skeleton
        if(structureName == 'skeleton'):
            ax1.plot(rotatedPosition[i,0,:],rotatedPosition[i,1,:],'bo')
            ax1.axis([-1.0,1.0,0.0,2.0])
            ax1.axis('equal')
    
    fig1.savefig(dirName+'RotatedDimersRe'+str(Re)+'.png')
    fig1.clf()
    
    return

def zipFiles(src,dst):
    zf = zipfile.ZipFile('%s.zip' % (dst), 'w', zipfile.ZIP_DEFLATED)
    abs_src = os.path.abspath(src)
    for dirname, subdirs, files in os.walk(src):
        for filename in files:
            absname = os.path.abspath(os.path.join(dirname, filename))
            arcname = absname[len(abs_src) + 1:]
            print('zipping %s as %s' % (os.path.join(dirname, filename),
                                        arcname))
            zf.write(absname, arcname)
    zf.close()

def MakeDirectory(directory,seed):
    if not os.path.exists(directory+'/'+str(seed)):
        os.makedirs(directory+'/'+str(seed))
    return
    
def main():
    
    #Number of swimmers
    Nbots = 2
    #Distance from COM
    #distFromCOM = 0.5
    #List of structure names
    structureNames = ['botup','botlow','skeleton']
    Reynolds = ['5','100']
    
    #Determine COM for structure
    COM = np.zeros((Nbots,2))
    
    #BOTUP
    #Read vertex file line by line and save the values in a list using the \n delimiter
    linesup = [line.strip() for line in open('botup'+str(1)+'.vertex')]
    #Break each list element into an array of numbers using the space delimiter
    linesup = [line.split() for line in linesup]
    nvertup = int(linesup[0][0])
    #BOTLOW
    #Read vertex file line by line and save the values in a list using the \n delimiter
    lineslow = [line.strip() for line in open('botlow'+str(1)+'.vertex')]
    #Break each list element into an array of numbers using the space delimiter
    lineslow = [line.split() for line in lineslow]
    nvertlow = int(lineslow[0][0])
    
    #Array of combined structure positions
    combinedPosition = np.zeros((2,nvertup+nvertlow))
    for i in range(1,nvertup+1):
        combinedPosition[0,i-1] = float(linesup[i][0])
        combinedPosition[1,i-1] = float(linesup[i][1])
    print('Xpos = ',combinedPosition[0,nvertup-1])
    print('Ypos = ',combinedPosition[1,nvertup-1])
    for i in range(nvertup+1,nvertup+nvertlow+1):
        combinedPosition[0,i-1] = float(lineslow[i-nvertup][0])
        combinedPosition[1,i-1] = float(lineslow[i-nvertup][1])
        if(i == nvertup+1):
            print('Xpos = ',combinedPosition[0,i-1])
            print('Ypos = ',combinedPosition[1,i-1])
    
    #Find the COM of our structure
    COM[:,0] = np.sum(combinedPosition[0,:])/float(nvertup+nvertlow)
    COMx = COM[0,0]
    print('COMx = ',COM[0,0])
    COM[:,1] = np.sum(combinedPosition[1,:])/float(nvertup+nvertlow)
    COMy = COM[0,1]
    print('COMy = ',COM[0,1])
    #return
    
    #Loops to create grid and make directories
    for Re in Reynolds:
        #Plot Dimer Locations
        fig0 = plt.figure(0)
        ax0 = fig0.add_subplot(1,1,1,aspect=1)
        ax0.set_title('Dimer Location Re'+str(Re))
        
        #Make Directory where vertex files and stuff will be stored
        MakeDirectory('../Structures','Re'+str(Re))
        dirName = '../Structures/Re'+str(Re)+'/'
        gridPosition, nvert, COM = DimerLocation(Nbots,structureNames,COM,dirName,Re,ax0)
        for i in range(len(structureNames)):
            RotateVertices(Nbots,nvert,structureNames[i],dirName,gridPosition,COM,Re)
        COM[0,0] = COMx
        COM[0,1] = COMy
        COM[1,0] = COMx
        COM[1,1] = COMy
        copyfile('2botSIBRe'+str(Re)+'.bsub',dirName+'2botSIBRe'+str(Re)+'.bsub')
        copyfile('input2dRe'+str(Re),dirName+'input2dRe'+str(Re))
        copyfile('main.C',dirName+'main.C')
        copyfile('Makefile',dirName+'Makefile')
        copyfile('update_springs.C',dirName+'update_springs.C')
        copyfile('update_springs.h',dirName+'update_springs.h')
        #Save PNG of dimer locations
        fig0.savefig(dirName+'DimerLocationRe'+str(Re)+'.png')
        fig0.clf()
        #Create Zip file of all .vertex, .spring, and .target files of interest
        #zipFiles(dirName,'Dimer'+str(i)+str(j))
            
    zipFiles('../Structures/','AYScattering')
    
#-------------------END MAIN---------------------
main()