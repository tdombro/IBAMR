#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 16:44:48 2018

@author: thomas
"""

Nbots = 9
Re = 100

pwd = "/Users/thomas/sfw/visitFiles/thomas/1bot/StandardIB/MultipleSwimmers/N256/Re5/Rotated/"

fluid = "localhost:"+pwd+"viz2D9botRe"+str(Re)+"/dumps.visit"
mesh = "localhost:"+pwd+"viz2D9botRe"+str(Re)+"/lag_data.visit"
dbs = (fluid,mesh)

OpenDatabase(dbs[0], 0)
DeleteAllPlots()

AddPlot("Pseudocolor", "Omega", 1, 1)
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
PseudocolorAtts.skewFactor = 1
PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
PseudocolorAtts.minFlag = 1
PseudocolorAtts.min = -200
PseudocolorAtts.maxFlag = 1
PseudocolorAtts.max = 200
PseudocolorAtts.centering = PseudocolorAtts.Nodal  # Natural, Nodal, Zonal
PseudocolorAtts.colorTableName = "difference"
PseudocolorAtts.invertColorTable = 0
PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
PseudocolorAtts.opacityVariable = ""
PseudocolorAtts.opacity = 1
PseudocolorAtts.opacityVarMin = 0
PseudocolorAtts.opacityVarMax = 1
PseudocolorAtts.opacityVarMinFlag = 0
PseudocolorAtts.opacityVarMaxFlag = 0
PseudocolorAtts.pointSize = 0.05
PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
PseudocolorAtts.pointSizeVarEnabled = 0
PseudocolorAtts.pointSizeVar = "default"
PseudocolorAtts.pointSizePixels = 2
PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
PseudocolorAtts.lineWidth = 0
PseudocolorAtts.tubeResolution = 10
PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.tubeRadiusAbsolute = 0.125
PseudocolorAtts.tubeRadiusBBox = 0.005
PseudocolorAtts.tubeRadiusVarEnabled = 0
PseudocolorAtts.tubeRadiusVar = ""
PseudocolorAtts.tubeRadiusVarRatio = 10
PseudocolorAtts.tailStyle = PseudocolorAtts.None  # None, Spheres, Cones
PseudocolorAtts.headStyle = PseudocolorAtts.None  # None, Spheres, Cones
PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.endPointRadiusAbsolute = 0.125
PseudocolorAtts.endPointRadiusBBox = 0.05
PseudocolorAtts.endPointResolution = 10
PseudocolorAtts.endPointRatio = 5
PseudocolorAtts.endPointRadiusVarEnabled = 0
PseudocolorAtts.endPointRadiusVar = ""
PseudocolorAtts.endPointRadiusVarRatio = 10
PseudocolorAtts.renderSurfaces = 1
PseudocolorAtts.renderWireframe = 0
PseudocolorAtts.renderPoints = 0
PseudocolorAtts.smoothingLevel = 0
PseudocolorAtts.legendFlag = 1
PseudocolorAtts.lightingFlag = 1
PseudocolorAtts.wireframeColor = (0, 0, 0, 0)
PseudocolorAtts.pointColor = (0, 0, 0, 0)
SetPlotOptions(PseudocolorAtts)

OpenDatabase(dbs[1], 0)
#SetDefaultAttributes
MeshAtts = MeshAttributes()
MeshAtts.legendFlag = 0
MeshAtts.lineStyle = MeshAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
MeshAtts.lineWidth = 0
MeshAtts.meshColor = (0, 0, 0, 255)
MeshAtts.meshColorSource = MeshAtts.Foreground  # Foreground, MeshCustom
MeshAtts.opaqueColorSource = MeshAtts.Background  # Background, OpaqueCustom
MeshAtts.opaqueMode = MeshAtts.Auto  # Auto, On, Off
MeshAtts.pointSize = 0.05
MeshAtts.opaqueColor = (255, 255, 255, 255)
MeshAtts.smoothingLevel = MeshAtts.None  # None, Fast, High
MeshAtts.pointSizeVarEnabled = 0
MeshAtts.pointSizeVar = "default"
MeshAtts.pointType = MeshAtts.Sphere  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
MeshAtts.showInternal = 0
MeshAtts.pointSizePixels = 5
MeshAtts.opacity = 1
SetDefaultPlotOptions(MeshAtts)

for i in range(1,Nbots+1):
    AddPlot("Mesh", "botlow"+str(i)+"_vertices", 1, 1)
    AddPlot("Mesh", "botup"+str(i)+"_vertices", 1, 1)
    
for i in range(1,Nbots+1):

    AddPlot("Mesh", "skeleton"+str(i)+"_mesh", 1, 1)
    MeshAtts = MeshAttributes()
    MeshAtts.legendFlag = 0
    MeshAtts.lineStyle = MeshAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    MeshAtts.lineWidth = 0
    MeshAtts.meshColor = (0, 255, 255, 255)
    MeshAtts.meshColorSource = MeshAtts.MeshCustom  # Foreground, MeshCustom
    SetPlotOptions(MeshAtts)

CreateDatabaseCorrelation("common",dbs,0)
DrawPlots()
w = GetWindowInformation()
print "Active time slider: %s" % w.timeSliders[w.activeTimeSlider]
SetActiveTimeSlider("common")
w = GetWindowInformation()
print "Active time slider: %s" % w.timeSliders[w.activeTimeSlider]

AnnotationAtts = AnnotationAttributes()
AnnotationAtts.axes2D.xAxis.title.visible = 0
AnnotationAtts.axes2D.yAxis.title.visible = 0
AnnotationAtts.userInfoFlag = 0
SetAnnotationAttributes(AnnotationAtts)

#Add Annotations
title = CreateAnnotationObject("Text2D")
title.text = "9bot seed=0 Re="+str(Re)
title.position = (0.35,0.95)
title.fontBold = 1
title.height = 0.03

ydist = CreateAnnotationObject("Text2D") 
ydist.text = "Y-dist (m)" 
ydist.position = (0.07, 0.92) 
ydist.height = 0.02
ydist.fontBold = 1  
print ydist


xdist = CreateAnnotationObject("Text2D") 
xdist.text = "X-dist (m)" 
xdist.position = (0.55, 0.09) 
xdist.height = 0.02
xdist.fontBold = 1  
print xdist

TimeSlider = CreateAnnotationObject("TimeSlider")
TimeSlider.position = (0.40,0.03)
TimeSlider.startColor = (255,0,0,255)
TimeSlider.endColor = (255,255,255,255)
TimeSlider.text = "Time=$time s"    

DrawPlots()

SaveSession(pwd+"visit9botRe"+str(Re)+".session")
