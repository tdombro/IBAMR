#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 14:16:38 2018

@author: thomas
"""

import numpy as np
import matplotlib.pyplot as plt

def main():
    
    nskeletons = 9
    
    #Plot Both Sets of Discs
    fig0 = plt.figure(0)
    ax0 = fig0.add_subplot(1,1,1,aspect=1)
    ax0.set_title('Non-Rotated Discs')
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(1,1,1,aspect=1)
    ax1.set_title('Rotated Discs')
    
    #Generate Random Angle!
    np.random.seed = 1123
    theta = np.zeros(nskeletons)
    for i in range(nskeletons):
        theta[i] = 2.0*np.pi*np.random.random()
        print('theta[%i] = %.3e' %(i,theta[i]))
    #return
    
    #Loop over all vertex files
    for i in range(1,nskeletons+1):
        #Read vertex file line by line and save the values in a list using the \n delimiter
        lines = [line.strip() for line in open('skeleton'+str(i)+'.vertex')]
        #Break each list element into an array of numbers using the space delimiter
        lines = [line.split() for line in lines]
        
        nvert = int(lines[0][0])
        print('nvert = ',nvert)
        
        #Allocate Arrays
        arrPosition = np.zeros((nskeletons,2,nvert))
        rotationMatrix = np.zeros((2,2))
        rotatedPosition = np.zeros((nskeletons,2,nvert))
        
        #Generate Random Angle!
        #np.random.seed = 1123
        #theta = 2.0*np.pi*np.random.rand(9)
        print('theta[%i] = %.3e'%(i,theta[i-1]))
        #Create rotation matrix
        rotationMatrix[0,0] = np.cos(theta[i-1])
        rotationMatrix[0,1] = -1.0*np.sin(theta[i-1])
        rotationMatrix[1,0] = np.sin(theta[i-1])
        rotationMatrix[1,1] = np.cos(theta[i-1])
        
        #Store vertex coord values
        for j in range(1,nvert+1):
            arrPosition[i-1,0,j-1] = lines[j][0]
            arrPosition[i-1,1,j-1] = lines[j][1]
        
        #Find center of mass
        #COM = 1/M * SUM(m_i * r_i)
        #M = nvert; m = 1
        #COM = 1/nvert * SUM(r_i)
        COMx = np.sum(arrPosition[i-1,0,:])/float(nvert)
        print('COMx = ',COMx)
        COMy = np.sum(arrPosition[i-1,1,:])/float(nvert)
        print('COMy = ',COMy)
        COM = np.array([COMx,COMy])
        
        #open newly rotated vertex file
        f = open('rotatedskeleton'+str(i)+'.vertex','w')
        f.write('%i\n'%nvert)
        
        #Rotate coord values about center of mass
        for j in range(nvert):
            rotatedPosition[i-1,:,j] = rotationMatrix.dot(arrPosition[i-1,:,j] - COM)
            rotatedPosition[i-1,:,j] += COM
            if(j == nvert - 1):
                f.write('%.5e %.5e'%(rotatedPosition[i-1,0,j],rotatedPosition[i-1,1,j]))
            else:
                f.write('%.5e %.5e\n' %(rotatedPosition[i-1,0,j],rotatedPosition[i-1,1,j]))
        f.close()    
        
        #Plot Skeleton
        ax0.plot(arrPosition[i-1,0,:],arrPosition[i-1,1,:],'ro')
        ax0.axis([-2.0,2.0,-4.0,4.0])
        ax0.axis('equal')
        
        #Plot Rotated Skeleton
        ax1.plot(rotatedPosition[i-1,0,:],rotatedPosition[i-1,1,:],'bo')
        ax1.axis([-2.0,2.0,-4.0,4.0])
        ax1.axis('equal')

    plt.show()
            
        
            
    
    
    
    
    
    
    
    
    
    
    
#-------------------END MAIN---------------------
main()