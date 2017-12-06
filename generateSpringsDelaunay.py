#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 15:23:16 2017

@author: thomas
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Ks = 1.0E5

def get_distance(vert1,vert2):
    distX = vert1[0] - vert2[0]
    distY = vert1[1] - vert2[1]
    dist = distX*distX + distY*distY
    return np.sqrt(dist)

def main():
    #Read radius file line by line and save the values in a list using the \n delimiter
    lines = [line.strip() for line in open('radius.txt')]
    #Break each list element into an array of numbers using the space delimiter
    lines = [line.split() for line in lines]

    #Determine # of vertices on surface of sphere
    totCount = 0
    index = np.zeros(0) #0 -> 236 values 1 -> 238 skip 7
    circleCount = 0

    for i in range(len(lines)):
        index = np.append(index,int(lines[i][0]))
        totCount += 1
    print('totCount = ',totCount)        
    
    #Array to store position (x,y,z) 2D and 3D
    #points = np.zeros((totCount,3)) #0 -> 236
    points = np.zeros((totCount,2)) #0 -> 236
    
    #Store position
    arrCount = 0
    for i in range(totCount):
        points[arrCount,0] = float(lines[i][1])
        points[arrCount,1] = float(lines[i][2])
        arrCount += 1
    print("arrCount = ",arrCount)
    print("circleCount = ",circleCount)
    
    #2D Delaunay Triangulation
    from scipy.spatial import Delaunay
    tri = Delaunay(points)
    
    #Plot Delaunay Triangulation
    plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
    plt.plot(points[:,0], points[:,1], 'o')
    plt.show()
    
    pairs = np.zeros((tri.simplices.size,2))
    pairCount = 0
    #Generate Pairs of vertices
    for simplex in tri.simplices:
        #Store Vertex pairs
        pairs[pairCount,:] = [simplex[0],simplex[1]]
        pairCount += 1
        
    print(pairCount)
    
    savedPairs = np.zeros((int(pairCount),2))
    boolPair = 1
    sPairCount = 0
    
    print(tri.simplices.size)
    
    #Report all non-repeating pairs
    for i in range(pairCount):
        for j in range(i+1, pairCount):
            if(j != i):
                if((pairs[i,0] == pairs[j,0] and pairs[i,1] == pairs[j,1]) 
                or (pairs[i,0] == pairs[j,1] and pairs[i,1] == pairs[j,0])):
                    boolPair = 0
                    break
            else:
                boolPair = 0
                break
        if(boolPair == 1):
            #convert index value used to find dist b/w 2 points
            savedPairs[sPairCount,:] = [pairs[i,0],pairs[i,1]]
            print("sPC[] = [%f, %f]\tsPC = %i" %(savedPairs[sPairCount,0], savedPairs[sPairCount,1],sPairCount))
            sPairCount += 1
        boolPair = 1
        
    sPx, sPy = savedPairs.shape
    print(sPx,sPy)
    
    #Generate spring file in a new file called disc_test.spring
    f = open('disc_test.spring','w')
    f.write('%i\n' % sPairCount)
    for i in range(sPairCount):
        #convert index value used to find vertex coordinates
        vert1 = [points[int(savedPairs[i,0]),0], points[int(savedPairs[i,0]),1]]
        vert2 = [points[int(savedPairs[i,1]),0], points[int(savedPairs[i,1]),1]]
        
        #Find distance b/w the 2 vertices
        dist = get_distance(vert1,vert2)
        print(dist)
        f.write('%i\t%i %f %f\n' % (int(savedPairs[i,0]),int(savedPairs[i,1]),Ks,dist))

#----------------_END MAIN_-------------------------------#
main()