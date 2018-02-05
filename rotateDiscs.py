#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 14:16:38 2018

@author: thomas
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import zipfile

#Plot Non-Rotated Discs
fig0 = plt.figure(0)
ax0 = fig0.add_subplot(1,1,1,aspect=1)
ax0.set_title('Non-Rotated Discs')

#Plot Rotated Discs
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(1,1,1,aspect=1)
ax1.set_title('Rotated Discs')

def GridSetup(Nbots,structureName,COM,distFromCOM):
    #Sets up a grid of swimmers to use in IBAMR simulation. Blocks them
    #grid starts at (0,0). Grid increments (0.9,0.9) for every additional row,column of swimmers
    #Set up to work for sqrt(Nbots) many swimmers
    
    #Index for bot #
    vertexNumber = 1
    
    #Determine # vertices in current structure
    #Read vertex file line by line and save the values in a list using the \n delimiter
    lines = [line.strip() for line in open(structureName+str(5)+'.vertex')]
    #Break each list element into an array of numbers using the space delimiter
    lines = [line.split() for line in lines]
    nvert = int(lines[0][0])
    print('nvert = ',nvert)
    
    #Allocate Arrays
    arrPosition = np.zeros((Nbots,2,nvert))
    
    for j in range(int(np.sqrt(Nbots))):
        for i in range(int(np.sqrt(Nbots))):
            #Read vertex file line by line and save the values in a list using the \n delimiter
            lines = [line.strip() for line in open(structureName+str(5)+'.vertex')]
            #Break each list element into an array of numbers using the space delimiter
            lines = [line.split() for line in lines]
            
            #Store vertex coord values
            f = open(structureName+'Grid'+str(vertexNumber)+'.vertex','w')
            f.write('%i\n'%nvert)
            for k in range(1,nvert+1):
                arrPosition[vertexNumber-1,0,k-1] = float(lines[k][0])
                arrPosition[vertexNumber-1,1,k-1] = float(lines[k][1])
                #Translate Vertex coordinates to fit in current grid cell
                arrPosition[vertexNumber-1,0,k-1] += distFromCOM*(1.0 + 2.0*i)
                arrPosition[vertexNumber-1,1,k-1] += distFromCOM*(1.0 + 2.0*j)
                if(k == nvert):
                    f.write('%.5e %.5e'%(arrPosition[vertexNumber-1,0,k-1],arrPosition[vertexNumber-1,1,k-1]))
                else:
                    f.write('%.5e %.5e\n'%(arrPosition[vertexNumber-1,0,k-1],arrPosition[vertexNumber-1,1,k-1]))
            #Translate COM
            if(structureName == 'skeleton'):
                COM[vertexNumber-1,0] += distFromCOM*(1.0 + 2.0*i)
                COM[vertexNumber-1,1] += distFromCOM*(1.0 + 2.0*j)
                print('COMx = ',COM[vertexNumber-1,0])
                print('COMy = ',COM[vertexNumber-1,1])

            f.close()

            #Plot Skeleton
            ax0.plot(arrPosition[vertexNumber-1,0,:],arrPosition[vertexNumber-1,1,:],'ro')
            ax0.axis([0.0,3.0,0.0,3.0])
            ax0.axis('equal')
                
            #Iterate vertex number
            vertexNumber += 1
    
    return (arrPosition,nvert,COM) 
    

def RotateVertices(Nbots,nvert,structureName,gridPosition,COM):  
    #Rotates the set up grid by a randomly generated angle
    
    #Allocate Arrays
    rotationMatrix = np.zeros((Nbots,2,2))
    theta = np.zeros(Nbots)
    rotatedPosition = np.zeros((Nbots,2,nvert))
    #COM = np.zeros((Nbots,2))
    
    #Generate Random Angles
    np.random.seed(0)
    for i in range(Nbots):
        theta[i] = 2.0*np.pi*np.random.rand(1)
        print('theta[%i] = %.3e' %(i,theta[i]))
        rotationMatrix[i,0,0] = np.cos(theta[i])
        rotationMatrix[i,0,1] = -1.0*np.sin(theta[i])
        rotationMatrix[i,1,0] = np.sin(theta[i])
        rotationMatrix[i,1,1] = np.cos(theta[i])     
        
    #Rotate Grid Structures around COM based off of above angle 
    for i in range(int(Nbots)):
        #open newly rotated vertex file
        f = open('RotatedDisc/'+structureName+str(i+1)+'.vertex','w')
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
        ax1.plot(rotatedPosition[i,0,:],rotatedPosition[i,1,:],'bo')
        ax1.axis([0.0,3.0,0.0,3.0])
        ax1.axis('equal')
    
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
    
def main():
    
    #Number of swimmers
    Nbots = 9
    #Distance from COM
    distFromCOM = 0.5
    #List of structure names
    structureNames = ['skeleton','botup','botlow']
    
    #Determine COM for structure
    COM = np.zeros((Nbots,2))
    
    #BOTUP
    #Read vertex file line by line and save the values in a list using the \n delimiter
    linesup = [line.strip() for line in open('botup'+str(5)+'.vertex')]
    #Break each list element into an array of numbers using the space delimiter
    linesup = [line.split() for line in linesup]
    nvertup = int(linesup[0][0])
    #BOTLOW
    #Read vertex file line by line and save the values in a list using the \n delimiter
    lineslow = [line.strip() for line in open('botlow'+str(5)+'.vertex')]
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
    print('COMx = ',COM[0,0])
    COM[:,1] = np.sum(combinedPosition[1,:])/float(nvertup+nvertlow)
    print('COMy = ',COM[0,1])
    #return
    
    for i in range(len(structureNames)):
        gridPosition, nvert, COM = GridSetup(Nbots,structureNames[i],COM,distFromCOM)
        RotateVertices(Nbots,nvert,structureNames[i],gridPosition,COM)
    
    ax0.plot(COM[:,0],COM[:,1],'ko')
    ax1.plot(COM[:,0],COM[:,1],'ko')
    for i in range(int(np.sqrt(Nbots)+1)):
        ax0.axvline(x=2.0*i*distFromCOM)
        ax0.axhline(y=2.0*i*distFromCOM)
        ax1.axvline(x=2.0*i*distFromCOM)
        ax1.axhline(y=2.0*i*distFromCOM)
        
    fig0.savefig('NonRotated.png')
    fig1.savefig('Rotated.png')
    #plt.show()
    
    zipFiles('RotatedDisc','RotatedDiscs')
    
#-------------------END MAIN---------------------
main()