#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri july 28 11:45:25 2017

@author: thomas
"""

import numpy as np

#Read vertex file line by line and save the values in a list using the \n delimiter
lines = [line.strip() for line in open('bot1low.vertex')]
#Break each list element into an array of numbers using the space delimiter
lines = [line.split() for line in lines]
print(lines[1])
print(len(lines))

posArr = np.zeros((3,len(lines)))
radius = np.zeros(int(len(lines)))
boundary = np.zeros(int(len(lines)))
#Find Position of each vertex
for i in range(1,len(lines)):
   posArr[0,i] = float(lines[i][0]) + 0.45
   posArr[1,i] = float(lines[i][1]) + 0.75
   posArr[2,i] = float(lines[i][2])
   #Determine Radius of each vertex
   radius[i] = np.sqrt(posArr[0,i]*posArr[0,i] + posArr[1,i]*posArr[1,i] + posArr[2,i]*posArr[2,i])
   #Is the vertex at the boundary?
   if(radius[i] + 1.0e-3 >= 0.15 ): 
       boundary[i] = 1
   else: 
       boundary[i] = 0
print(posArr[0,1])
print(radius)
print(boundary)
#Record Radius in a new file called radius.text
f = open('radius.txt','w')
for i in range(1,np.size(radius)):
    f.write("%i %f %f %f %f %i\n"% (i-1,posArr[0,i],posArr[1,i],posArr[2,i],radius[i],boundary[i]))
f.close()