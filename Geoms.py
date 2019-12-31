# -*- coding: utf-8 -*-
#For 3D molecule visualization using vector mathematics

from vpython import *
from copy import *
import math
from tkinter import *
import re


def computerOrthogonals(vector):
    #There are two possible perpendicular vectors, and we will find both by 
    #   switching and negating the original vector's x + y components. The original
    #   vector is <x, y>. The left will be <y, -x> and the right will be <-y, x>
    leftPerpVector, rightPerpVector = copy(vector), copy(vector)
    #computing the left vector
    leftPerpVector[0], leftPerpVector[1] = leftPerpVector[1], leftPerpVector[0]
    leftPerpVector[1] = -leftPerpVector[1]
    rightPerpVector[0], rightPerpVector[1] = rightPerpVector[1], rightPerpVector[0]
    rightPerpVector[0] = -rightPerpVector[0]
    return leftPerpVector, rightPerpVector

#Refer to the color code chart
def figureColor(element):
    element = re.split('(\d+)', element)[0]
    if element == 'H':
        return vector(1, 1, 1) #White
    elif element == 'C':
        return vector(0.2, 0.2, 0.2) #Black
    elif element == 'F':
        return vector(1, 0, 1) #Magenta
    elif element == 'Cl':
        return vector(1, 0.8, 1) #Green
    elif element == 'Br':
        return vector(0.6, 0.298, 0) #Brown
    elif element == 'I':
        return vector(1, 0.6, 0) #Orange
    elif element == 'O':
        return vector(1, 0, 0) #Red
    elif element == 'S':
        return vector(1, 1, 0) #Yellow
    elif element == 'Se':
        return vector(0.4, 0.4, 0.4) #Gray
    elif element == 'Te':
        return vector(1, 0.7, 0.2) #Coppery
    elif element == 'N':
        return vector(0, 0, 1) #blue
    elif element == 'P':
        return vector(1, 0.8, 0.898)
    elif element == 'As':
        return vector(0.6, 1, 1)
    elif element == 'Sb':
        return vector(0, 0.4, 0)
    elif element == 'Si':
        return vector(0.8, 1, 0.8)
    elif element == 'Ge':
        return vector(1, 0.5, 0.937)
    elif element == 'B':
        return vector(0.1529, 0.412, 0.439)
    elif element == 'Xe':
        return vector(0.76, 0.72, 0.455)
    elif element == 'Kr':
        return vector(0.098, 0.475, 0.00392)


#Inverts a vector
def invert(vectorObj):
    x, y, z = vectorObj.x, vectorObj.y, vectorObj.z
    x, y, z = -x, -y, -z
    return vector(x, y, z)

def captureScreen(filename):
    file = filename + '.png'
    scene.capture(file)
    
#Clears the screen b/w updates of the view
def clearScreen(filename, weight):
    scene.width = 1900
    scene.height = 950
    scene.caption = '3D Molecule: %s    Molecule Weight: %0.2f g/mol    Right-Click + drag to rotate, mouse wheel to zoom.    ' % (filename, weight)
    button(bind = lambda: captureScreen(filename), text = 'Screenshot!')
    scene.append_to_caption('\t\t')
    for obj in scene.objects:
        obj.visible = False
        del obj

#Finds the required objects from the scene
def findNecessaryObjects(currAtom, prevAtom):
    currAtomObj, prevAtomObj = None, None
    for obj in scene.objects:
        if isinstance(obj, sphere):
            if obj.label == currAtom: currAtomObj = obj
            elif obj.label == prevAtom: prevAtomObj = obj
    return currAtomObj, prevAtomObj

#SurroundingAtoms is a list of vector objects
def drawBonds(centAtomPos, surroundingAtoms, bondLst):
    for i in range(len(surroundingAtoms)):
        difference = surroundingAtoms[i] - centAtomPos
        if float(difference.z) == 0: #In the plane somewhere
            if bondLst[i] == 'tb':
                listVector = [difference.x, difference.y]
                leftPerpVector, rightPerpVector = computerOrthogonals(listVector)
                leftPerpVec = vector(leftPerpVector[0], leftPerpVector[1], 0)
                rightPerpVec = vector(rightPerpVector[0], rightPerpVector[1], 0)
                leftPerpNorm, rightPerpNorm = leftPerpVec.norm(), rightPerpVec.norm()
                leftPerpNormAdded, rightPerpNormAdded = (0.5) * leftPerpNorm, (0.5) * rightPerpNorm
                startPos1 = vector(centAtomPos.x + leftPerpNormAdded.x, centAtomPos.y + leftPerpNormAdded.y, centAtomPos.z + leftPerpNormAdded.z)
                startPos2 = vector(centAtomPos.x + rightPerpNormAdded.x, centAtomPos.y + rightPerpNormAdded.y, centAtomPos.z + rightPerpNormAdded.z)
                startPosOrig = centAtomPos
                bond1 = cylinder(pos = startPos1, axis = difference, radius = 0.15)
                bond2 = cylinder(pos = startPos2, axis = difference, radius = 0.15)
                bond3 = cylinder(pos = startPosOrig, axis = difference, radius = 0.15)
            elif bondLst[i] == 'db':
                listVector = [difference.x, difference.y]
                leftPerpVector, rightPerpVector = computerOrthogonals(listVector)
                leftPerpVec = vector(leftPerpVector[0], leftPerpVector[1], 0)
                rightPerpVec = vector(rightPerpVector[0], rightPerpVector[1], 0)
                leftPerpNorm, rightPerpNorm = leftPerpVec.norm(), rightPerpVec.norm()
                leftPerpNormAdded, rightPerpNormAdded = (0.5) * leftPerpNorm, (0.5) * rightPerpNorm
                startPos1 = vector(centAtomPos.x + leftPerpNormAdded.x, centAtomPos.y + leftPerpNormAdded.y, centAtomPos.z + leftPerpNormAdded.z)
                startPos2 = vector(centAtomPos.x + rightPerpNormAdded.x, centAtomPos.y + rightPerpNormAdded.y, centAtomPos.z + rightPerpNormAdded.z)
                bond1 = cylinder(pos = startPos1, axis = difference, radius = 0.20)
                bond2 = cylinder(pos = startPos2, axis = difference, radius = 0.20)
            elif bondLst[i] == 'sb':
                bond = cylinder(pos = centAtomPos, axis = difference, radius = 0.25)
        else: #We're somewhere in 3D space
            x2, y2 = 1, 0 #Arbitrary values for x2, y2
            x3, y3 = -1, 0
            x1, y1, z1 = difference.x, difference.y, difference.z
            z2 = (-x1 * x2 - y1 * y2) / z1
            z3 = (-x1 * x3 - y1 * y3) / z1
            perpVector1, perpVector2 = vector(x2, y2, z2), vector(x3, y3, z3)
            perpVector1Norm, perpVector2Norm = perpVector1.norm(), perpVector2.norm()
            perpVec1Added = (0.5) * perpVector1Norm
            perpVec2Added = (0.5) * perpVector2Norm
            if bondLst[i] == 'tb':
                startPos1 = vector(centAtomPos.x + perpVec1Added.x, centAtomPos.y + perpVec1Added.y, centAtomPos.z + perpVec1Added.z)
                startPos2 = vector(centAtomPos.x + perpVec2Added.x, centAtomPos.y + perpVec2Added.y, centAtomPos.z + perpVec2Added.z)
                startPosOrig = centAtomPos
                bond1 = cylinder(pos = startPos1, axis = difference, radius = 0.15)
                bond2 = cylinder(pos = startPos2, axis = difference, radius = 0.15)
                bond3 = cylinder(pos = startPosOrig, axis = difference, radius = 0.15)
            elif bondLst[i] == 'db':
                startPos1 = vector(centAtomPos.x + perpVec1Added.x, centAtomPos.y + perpVec1Added.y, centAtomPos.z + perpVec1Added.z)
                startPos2 = vector(centAtomPos.x + perpVec2Added.x, centAtomPos.y + perpVec2Added.y, centAtomPos.z + perpVec2Added.z)
                bond1 = cylinder(pos = startPos1, axis = difference, radius = 0.20)
                bond2 = cylinder(pos = startPos2, axis = difference, radius = 0.20)
            elif bondLst[i] == 'sb':
                bond = cylinder(pos = centAtomPos, axis = difference, radius = 0.25)


#SN = 4
def startingTetrahedral(dist, atomList, bondLst):
    centerAtom = sphere(color = figureColor(atomList[0]))
    centerAtom.label = atomList[0] #We may not have to do this copy pasting later on, but that can be changed later
    atomListCopy = copy(atomList)
    atomListCopy = atomListCopy[1:]
    x1, x2 = (0.5) * dist, (-0.5) * dist
    y1, y2 = (0.5) * dist, (-0.5) * dist
    z1, z2 = (0.5) * dist, (-0.5) * dist
    possiblePoints = [vector(x2, y2, z2), vector(x1, y1, z2), vector(x2, y1, z1), vector(x1, y2, z1)]
    place = 0
    for atoms in atomListCopy:
        point = possiblePoints[place]
        col = figureColor(atoms)
        newSphere = sphere(pos = point, color = col)
        newSphere.label = atoms
        place += 1
    takenPossiblePoints = possiblePoints[:place]
    drawBonds(vector(0, 0, 0), takenPossiblePoints, bondLst)

#Current Atom is the one you're considering, previous atom is the one that's already on the screen
def addingTetrahedral(currAtomLabel, prevAtomLabel, arbitraryDist, newConnectingAtoms, bondLst):
    diagonalDistance = (2 ** 0.5) * arbitraryDist
    midPointDistance = ((2 ** 0.5) / 2) * arbitraryDist
    angle = math.radians(109.5)
    angle2 = math.radians(120)
    newVectors = []
    currAtomObj, prevAtomObj = findNecessaryObjects(currAtomLabel, prevAtomLabel)
    startVec, endVec = prevAtomObj.pos, currAtomObj.pos
    resultVector = endVec - startVec #Calculate the 3D vector b/w the original and new atom
    if float(resultVector.z) == 0: #We are in the xy plane
        result2DVec = resultVector.norm()
        lstVector = [result2DVec.x, result2DVec.y]
        leftOrtho, rightOrtho = computerOrthogonals(lstVector)
        orthogonalVec = vector(rightOrtho[0], rightOrtho[1], 0)
        firstDistal = rotate(resultVector, angle, orthogonalVec)
        secondDistal, thirdDistal = rotate(firstDistal, angle2, resultVector), rotate(firstDistal, 2 * angle2, resultVector)
        firstDistal, secondDistal, thirdDistal = invert(firstDistal), invert(secondDistal), invert(thirdDistal)
        x1, y1, z1 = endVec.x + firstDistal.x, endVec.y + firstDistal.y, endVec.z + firstDistal.z
        x2, y2, z2 = endVec.x + secondDistal.x, endVec.y + secondDistal.y, endVec.z + secondDistal.z
        x3, y3, z3 = endVec.x + thirdDistal.x, endVec.y + thirdDistal.y, endVec.z + thirdDistal.z
        newVectors.append(vector(x1, y1, z1)), newVectors.append(vector(x2, y2, z2)), newVectors.append(vector(x3, y3, z3))
        place = 0
        for atoms in newConnectingAtoms:
            pointVec = newVectors[place]
            color = figureColor(atoms)
            newAtom = sphere(pos = pointVec, color = color)
            newAtom.label = atoms
            place += 1
        takenPossiblePoints = newVectors[:place]
        drawBonds(endVec, takenPossiblePoints, bondLst)
        #Do this using more vector calculus
    elif float(resultVector.x) == float(resultVector.y) == 0: #We are along the z axis
        axis = vector(0, 1, 0) #Just rotate around the y axis
        firstDistal = rotate(resultVector, angle, axis)
        secondDistal, thirdDistal = rotate(firstDistal, angle2, resultVector), rotate(firstDistal, 2 * angle2, resultVector)
        firstDistal, secondDistal, thirdDistal = invert(firstDistal), invert(secondDistal), invert(thirdDistal)
        x1, y1, z1 = endVec.x + firstDistal.x, endVec.y + firstDistal.y, endVec.z + firstDistal.z
        x2, y2, z2 = endVec.x + secondDistal.x, endVec.y + secondDistal.y, endVec.z + secondDistal.z
        x3, y3, z3 = endVec.x + thirdDistal.x, endVec.y + thirdDistal.y, endVec.z + thirdDistal.z
        newVectors.append(vector(x1, y1, z1)), newVectors.append(vector(x2, y2, z2)), newVectors.append(vector(x3, y3, z3))
        place = 0
        for atoms in newConnectingAtoms:
            pointVec = newVectors[place]
            color = figureColor(atoms)
            newAtom = sphere(pos = pointVec, color = color)
            newAtom.label = atoms
            place += 1
        takenPossiblePoints = newVectors[:place]
        drawBonds(endVec, takenPossiblePoints, bondLst)
    else: #Somewhere in 3D space
        result2DVec = vector(resultVector.x, resultVector.y, 0)
        normal2DVec = result2DVec.norm() #Generates the unit vector
        vecToBeAdded = diagonalDistance * normal2DVec
        endX, endY = startVec.x + vecToBeAdded.x, startVec.y + vecToBeAdded.y #We found the directly opposite vertex
        newVectors.append(vector(endX, endY, startVec.z)) #Found one of the good points
        #Now calculating the other two vertices
        midVec = midPointDistance * normal2DVec
        midX, midY = startVec.x + midVec.x, startVec.y + midVec.y
        tempVec = [midX, midY]
        rightPerpVec, leftPerpVec = computerOrthogonals(tempVec) #Two vectors of the form [x, y]
        rightPerpVector, leftPerpVector = vector(rightPerpVec[0], rightPerpVec[1], 0), vector(leftPerpVec[0], leftPerpVec[1], 0)
        rightPerpVecNorm, leftPerpVecNorm = rightPerpVector.norm(), leftPerpVector.norm()
        vert1X, vert1Y = midX + (midPointDistance * rightPerpVecNorm.x), midY + (midPointDistance * rightPerpVecNorm.y)
        vert2X, vert2Y = midX + (midPointDistance * leftPerpVecNorm.x), midY + (midPointDistance * leftPerpVecNorm.y)
        vert1Z, vert2Z = None, None
        if startVec.z < endVec.z:
            vert1Z = vert2Z = startVec.z + arbitraryDist
        elif startVec.z > endVec.z:
            vert1Z = vert2Z = startVec.z - arbitraryDist
        newVectors.append(vector(vert1X, vert1Y, vert1Z)), newVectors.append(vector(vert2X, vert2Y, vert2Z))
        place = 0
        for atoms in newConnectingAtoms:
            pointVec = newVectors[place]
            color = figureColor(atoms)
            newAtom = sphere(pos = pointVec, color = color)
            newAtom.label = atoms
            place += 1
        takenPossiblePoints = newVectors[:place]
        drawBonds(endVec, takenPossiblePoints, bondLst)

#SN = 3
def startingTrigPlanar(dist, atomList, bondLst):
    centerAtom = sphere(color = figureColor(atomList[0]))
    centerAtom.label = atomList[0]
    atomListCopy = copy(atomList)
    atomListCopy = atomListCopy[1:]
    x1, y1, z1 = (0.5) * dist, 0, 0 #The first atom
    firstVec = vector(x1, y1, z1)
    possiblePoints = [firstVec]
    axis = vector(0, 0, 1) #Rotates around the z axis
    angle = math.radians(120)
    secondVec = rotate(firstVec, angle, axis)
    thirdVec = rotate(firstVec, 2 * angle, axis)
    possiblePoints.append(secondVec), possiblePoints.append(thirdVec)
    place = 0
    for atoms in atomListCopy:
        point = possiblePoints[place]
        color = figureColor(atoms)
        newAtom = sphere(pos = point, color = color)
        newAtom.label = atoms
        place += 1
    takenPossiblePoints = possiblePoints[:place]
    drawBonds(vector(0, 0, 0), takenPossiblePoints, bondLst)

#Easy way to do this would be to use vector rotation
def addingTrigPlanar(currAtomLabel, prevAtomLabel, arbitraryDist, newConnectingAtoms, bondLst):
    length = (0.5) * arbitraryDist
    newVectors = []
    currAtomObj, prevAtomObj = findNecessaryObjects(currAtomLabel, prevAtomLabel)
    startVec, endVec = prevAtomObj.pos, currAtomObj.pos
    difference = endVec - startVec 
    axis = None
    if float(difference.z) == 0: #We are in the plane
        axis = vector(0, 0, 1) #Rotate about the z axis
    else: #Somewhere in 3D space, not in the plane
        x2, y2 = 1, 0 #Arbitrary values for x2, y2
        x1, y1, z1 = difference.x, difference.y, difference.z
        z2 = (-x1 * x2 - y1 * y2) / z1
        axis = vector(x2, y2, z2)
    angle = math.radians(120)
    firstVec, secondVec = rotate(difference, angle, axis), rotate(difference, 2 * angle, axis)
    firstVec, secondVec = invert(firstVec),invert(secondVec)
    firstVecNorm, secondVecNorm = firstVec.norm(), secondVec.norm()
    firstX, firstY, firstZ = endVec.x + (length * firstVecNorm.x), endVec.y + (length * firstVecNorm.y), endVec.z + (length * firstVecNorm.z)
    secX, secY, secZ = endVec.x + (length * secondVecNorm.x), endVec.y + (length * secondVecNorm.y), endVec.z + (length * secondVecNorm.z)
    newVectors.append(vector(firstX, firstY, firstZ)), newVectors.append(vector(secX, secY, secZ))
    place = 0
    for atoms in newConnectingAtoms:
        pointVec = newVectors[place]
        color = figureColor(atoms)
        newAtom = sphere(pos = pointVec, color = color)
        newAtom.label = atoms
        place += 1
    takenPossiblePoints = newVectors[:place]
    drawBonds(endVec, takenPossiblePoints, bondLst)

#SN = 5
def startingTrigBiPyramid(dist, atomList, bondLst):
    length = (0.5) * dist
    centerAtom = sphere(color = figureColor(atomList[0]))
    centerAtom.label = atomList[0]
    pointOne, pointTwo = vector(0, 0, length), vector(0, 0, -length) #The two axial atoms of the trigonal bipyramid
    pointThree = vector(length, 0, 0) #The first planar atom
    possiblePoints = [pointOne, pointTwo, pointThree]
    angle = math.radians(120)
    axis = vector(0, 0, 1) #The z axis
    pointFour, pointFive = rotate(pointThree, angle, axis), rotate(pointThree, 2 * angle, axis)
    possiblePoints.append(pointFour), possiblePoints.append(pointFive)
    place = 0
    for atom in atomList[1:]:
        pointVec = possiblePoints[place]
        color = figureColor(atom)
        newAtom = sphere(pos = pointVec, color = color)
        newAtom.label = atom
        place += 1
    takenPossiblePoints = possiblePoints[:place]
    drawBonds(vector(0, 0, 0), takenPossiblePoints, bondLst)
        
def addingTrigBiPyramid(currAtomLabel, prevAtomLabel, arbitraryDist, newConnectingAtoms, bondLst):
    length = arbitraryDist * 0.5
    newVectors = [] #Expected length of four
    currAtomObj, prevAtomObj = findNecessaryObjects(currAtomLabel, prevAtomLabel)
    startVec, endVec = prevAtomObj.pos, currAtomObj.pos
    difference = endVec - startVec
    unitDifference = difference.norm() #The unit vector in the direction of the bond
    newSegment = unitDifference * length
    secX, secY, secZ = endVec.x + newSegment.x, endVec.y + newSegment.y, endVec.z + newSegment.z
    newVectors.append(vector(secX, secY, secZ)) #The directly opposite distal atom
    axis = difference #We are rotating about the central axis
    firstEqVec = None #The first equatorial atom
    if float(difference.z) == 0: #In the plane somewhere
        firstEqVec = vector(0, 0, 1) #Essentially the z-axis
    else:
        x2, y2 = 1, 0
        x1, y1, z1 = difference.x, difference.y, difference.z
        z2 = (-x1 * x2 - y1 * y2) / z1
        firstEqVec = vector(x2, y2, z2)
    angle = math.radians(120)
    secEqVec, thrEqVec = rotate(firstEqVec, angle, axis), rotate(firstEqVec, 2 * angle, axis) #The other two equatorial atoms
    firstEqNorm, secEqNorm, thrEqNorm = firstEqVec.norm(), secEqVec.norm(), thrEqVec.norm()
    firstNormLen, secNormLen, thrNormLen = firstEqNorm * length, secEqNorm * length, thrEqNorm * length
    for vecs in [firstNormLen, secNormLen, thrNormLen]:
        x, y, z = endVec.x + vecs.x, endVec.y + vecs.y, endVec.z + vecs.z
        newVectors.append(vector(x, y, z))
    place = 0
    for atom in newConnectingAtoms:
        pointVec = newVectors[place]
        color = figureColor(atom)
        newAtom = sphere(pos = pointVec, color = color)
        newAtom.label = atom
        place += 1
    takenPossiblePoints = newVectors[:place]
    drawBonds(endVec, takenPossiblePoints, bondLst)
    
#SN = 2
def startingLinear(dist, atomList, bondLst):
    length = (0.5) * dist
    centerAtom = sphere(color = figureColor(atomList[0]))
    centerAtom.label = atomList[0]
    atomListCopy = copy(atomList)
    atomListCopy = atomListCopy[1:]
    x1, y1, z1 = length, 0, 0
    x2, y2, z2 = -length, 0, 0
    possiblePoints = [vector(x1, y1, z1), vector(x2, y2, z2)]
    place = 0
    for atoms in atomListCopy:
        point = possiblePoints[place]
        color = figureColor(atoms)
        newAtom = sphere(pos = point, color = color)
        newAtom.label = atoms
        place += 1
    takenPossiblePoints = possiblePoints[:place]
    drawBonds(vector(0, 0, 0), takenPossiblePoints, bondLst)

def addingLinear(currAtomLabel, prevAtomLabel, arbitraryDist, newConnectingAtoms, bondLst):
    length = (0.5) * arbitraryDist
    currAtomObj, prevAtomObj = findNecessaryObjects(currAtomLabel, prevAtomLabel)
    startVec, endVec = prevAtomObj.pos, currAtomObj.pos
    difference = endVec - startVec
    difference = difference.norm()
    x2, y2, z2 = endVec.x + (length * difference.x), endVec.y + (length * difference.y), endVec.z + (length * difference.z)
    newVector = vector(x2, y2, z2)
    color = figureColor(newConnectingAtoms[0])
    newAtom = sphere(pos = newVector, color = color)
    newAtom.label = newConnectingAtoms[0] #There should only be one atom left in linear
    drawBonds(endVec, [newVector], bondLst)

#TODO: Make it so that if you have four atoms, you start in a planar formation
#SN = 6
def startingOctahedral(dist, atomList, bondLst):
    length = (0.5) * dist
    color = figureColor(atomList[0])
    centerAtom = sphere(color = color)
    centerAtom.label = atomList[0]
    possiblePoints = [vector(0, 0, -length), vector(0, 0, length)] #Initiate with the axial atoms first
    firstEq = vector(length, 0, 0) #The first equatorial atom
    possiblePoints.append(firstEq)
    angle, axis = math.radians(90), vector(0, 0, 1)
    secEq, thrEq, frthEq = rotate(firstEq, 2 * angle, axis), rotate(firstEq, angle, axis), rotate(firstEq, 3 * angle, axis)
    possiblePoints.append(secEq), possiblePoints.append(thrEq), possiblePoints.append(frthEq)
    place = 0
    for atom in atomList[1:]:
        newPos = possiblePoints[place]
        color = figureColor(atom)
        newAtom = sphere(pos = newPos, color = color)
        newAtom.label = atom
        place += 1
    takenPossiblePoints = possiblePoints[:place]
    drawBonds(vector(0, 0, 0), takenPossiblePoints, bondLst)

def addingOctahedral(currAtomLabel, prevAtomLabel, arbitraryDist, newConnectingAtoms, bondLst):
    length = (0.5) * arbitraryDist
    newVectors = []
    currAtomObj, prevAtomObj = findNecessaryObjects(currAtomLabel, prevAtomLabel)
    startVec, endVec = prevAtomObj.pos, currAtomObj.pos
    difference = endVec - startVec
    normedDifference = difference.norm()
    x1, y1, z1 = endVec.x + (length * normedDifference.x), endVec.y + (length * normedDifference.y), endVec.z + (length * normedDifference.z)
    newVectors.append(vector(x1, y1, z1))
    axis = difference
    firstEqVec = None
    if float(difference.z) == 0: #We are in the plane
        firstEqVec = vector(0, 0, 1)
    else:
        x2, y2 = 1, 0
        x1, y1, z1 = difference.x, difference.y, difference.z
        z2 = (-x1 * x2 - y1 * y2) / z1
        firstEqVec = vector(x2, y2, z2)
    angle = math.radians(90)
    secEq, thrEq, frthEq = rotate(firstEqVec, 2 * angle, axis), rotate(firstEqVec, angle, axis), rotate(firstEqVec, 3 * angle, axis)
    firstEqNorm, secEqNorm, thrEqNorm, frthEqNorm = firstEqVec.norm(), secEq.norm(), thrEq.norm(), frthEq.norm()
    firstNormLen, secNormLen, thrNormLen, frthNormLen = firstEqNorm * length, secEqNorm * length, thrEqNorm * length, frthEqNorm * length
    for vecs in [firstNormLen, secNormLen, thrNormLen, frthNormLen]:
        x, y, z = endVec.x + vecs.x, endVec.y + vecs.y, endVec.z + vecs.z
        newVectors.append(vector(x, y, z))
    place = 0
    for atom in newConnectingAtoms:        
        pointVec = newVectors[place]
        color = figureColor(atom)
        newAtom = sphere(pos = pointVec, color = color)
        newAtom.label = atom
        place += 1
    takenPossiblePoints = newVectors[:place]
    drawBonds(endVec, takenPossiblePoints, bondLst)

#Dealing with None type geometries
def placingNone(fullAtomLst, dist, bondLst):
    possibleVecs = [vector(0, 0, 0), vector((0.5) * dist, 0, 0)]
    place = 0
    for atom in fullAtomLst:
        pointVec = possibleVecs[place]
        color = figureColor(atom)
        newAtom = sphere(pos = pointVec, color = color)
        newAtom.label = atom
        place += 1
    drawBonds(vector(0, 0, 0), possibleVecs[1:], bondLst)
        
        
#Unpacks the tuples
def unpackAtomBondList(lst):
    atomLst, bondLst = [], []
    for (atom, bond) in lst:
        atomLst.append(atom), bondLst.append(bond)
    return atomLst, bondLst

#Code for drawing rings
#The length of atom lst should be equivalent to the length of the possible points list
def startingRings(radius, atomLst, matrix):
    numVertices = len(atomLst)
    possiblePoints = []
    sectAngle = (2 * math.pi) / numVertices
    for i in range(numVertices):
        angle = i * sectAngle
        newX, newY = (radius * math.cos(angle)), (radius * math.sin(angle))
        possiblePoints.append(vector(newX, newY, 0))
    for j in range(len(atomLst)):
        atomID, point = atomLst[j], possiblePoints[j]
        color = figureColor(atomID)
        newAtom = sphere(pos = point, color = color)
        newAtom.label = atomID
    for i in range(0, len(atomLst), 2): #Step by twos
        atom = atomLst[i]
        connectingList = matrix[atom]
        stringIndex = None
        for i in range(len(matrix[atom])):
            if type(matrix[atom][i]) == str and 'ring' in matrix[atom][i]:
                stringIndex = i
        atomsBonded = connectingList[:stringIndex]
        atomList, bondList = unpackAtomBondList(atomsBonded)
        firstAtom, secondAtom = atomList[0], atomList[1]
        firstAtomObj, secondAtomObj = findNecessaryObjects(firstAtom, secondAtom)
        currAtomObj = None
        for obj in scene.objects:
            if isinstance(obj, sphere):
                if obj.label == atom: currAtomObj = obj
        difference1 = firstAtomObj.pos - currAtomObj.pos
        difference2 = secondAtomObj.pos - currAtomObj.pos
        bond1 = cylinder(pos = currAtomObj.pos, axis = difference1, radius = 0.25)
        bond2 = cylinder(pos = currAtomObj.pos, axis = difference2, radius = 0.25)
        
        
def addingRings(currAtom, prevAtom, ringAtoms, radius):
    currAtomObj, prevAtomObj = None, None
    for obj in scene.objects:
        if isinstance(obj, sphere):    
            if obj.label == currAtom: currAtomObj = obj
            elif obj.label == prevAtom: prevAtomObj = obj
    endVec, startVec = currAtomObj.pos, prevAtomObj.pos
    difference = endVec - startVec
    diffNorm = difference.norm()
    distanceToCenter = radius * diffNorm
    #Coordinates of the center
    cX, cY, cZ = (endVec.x + distanceToCenter.x, endVec.y + distanceToCenter.y, endVec.z + distanceToCenter.z)
    numVertices = len(ringAtoms)
    sectorAngle = (2 * math.pi) / numVertices
    axis = None
    if float(difference.z) == 0: #in the xy plane somewhere
        axis = vector(0, 0, 1)
    else:
        x2, y2 = 1, 0
        x1, y1, z1 = difference.x, difference.y, difference.z
        z2 = (-x1 * x2 - y1 * y2) / z1
        axis = vector(x2, y2, z2)
    newVectors = []
    for i in range(numVertices):
        rotVec = rotate(diffNorm, i * sectorAngle, axis)
        rotVec = invert(rotVec)
        rotVecToBeAdded = (radius * rotVec)
        newX, newY, newZ = cX + rotVecToBeAdded.x, cY + rotVecToBeAdded.y, cZ + rotVecToBeAdded.z
        newVectors.append(vector(newX, newY, newZ))
    totalVectors = copy(newVectors)
    newVectors.pop(0) #No need to redraw the first atom
    atomPlace = ringAtoms.index(currAtom)
    atomPlace += 1
    #Algorithm builds COUNTER CLOCKWISE!! keep that in mind when fixing the order of atom labelling
    for vec in newVectors:
        color = figureColor(ringAtoms[atomPlace % len(ringAtoms)])
        newAtom = sphere(pos = vec, color = color)
        newAtom.label = ringAtoms[atomPlace % len(ringAtoms)]
        atomPlace += 1
    for i in range(len(totalVectors)):
        startPos, endPos = totalVectors[i], totalVectors[i - 1]
        bond = cylinder(pos = startPos, axis = endPos - startPos, radius = 0.25)
        

#This works...for now
def addingAtomsToRingsOne(adjacentAtoms, currAtom, atomToBeAdded, dist, bondType):
    dist = dist / 1.5
    leftAtom, rightAtom = adjacentAtoms[0], adjacentAtoms[1]
    currAtomVec, leftAtomVec, rightAtomVec = None, None, None
    for obj in scene.objects:
        if isinstance(obj, sphere):
            if obj.label == leftAtom: leftAtomVec = obj.pos
            elif obj.label == rightAtom: rightAtomVec = obj.pos
            elif obj.label == currAtom: currAtomVec = obj.pos
    leftBranch = leftAtomVec - currAtomVec
    rightBranch = rightAtomVec - currAtomVec
    dotProduct = dot(leftBranch, rightBranch)
    #A dot B = mag(A) * mag(B) * cos(theta)
    cosTheta = dotProduct / (mag(leftBranch) * mag(rightBranch))
    angle = math.acos(cosTheta) #Angle should be in radians
    reflexAngle = (2 * math.pi) - angle
    halfReflex = (0.5) * reflexAngle
    axis = cross(leftBranch, rightBranch)
    tempPoint = rotate(leftBranch, halfReflex, -axis)
    normTempVector = tempPoint.norm()
    tempVecToBeAdded = (dist / 2) * normTempVector
    newVec = vector(currAtomVec.x + tempVecToBeAdded.x, currAtomVec.y + tempVecToBeAdded.y, currAtomVec.z + tempVecToBeAdded.z)
    color = figureColor(atomToBeAdded[0])
    newAtom = sphere(pos = newVec, color = color)
    newAtom.label = atomToBeAdded[0]
    drawBonds(currAtomVec, [newVec], bondType)

def addingAtomsToRingsTwo(adjacentAtoms, currAtom, atomToBeAdded, dist, bondLst):
    dist = dist / 3
    leftAtom, rightAtom = adjacentAtoms[0], adjacentAtoms[1]
    currAtomVec, leftAtomVec, rightAtomVec = None, None, None
    for obj in scene.objects:
        if isinstance(obj, sphere):          
            if obj.label == leftAtom: leftAtomVec = obj.pos
            elif obj.label == rightAtom: rightAtomVec = obj.pos
            elif obj.label == currAtom: currAtomVec = obj.pos        
    leftBranch = leftAtomVec - currAtomVec
    rightBranch = rightAtomVec - currAtomVec
    perpVector = cross(leftBranch, rightBranch)
    perpNorm = perpVector.norm()
    perpVecToBeAdded = (dist / 2) * perpNorm
    dotProduct = dot(leftBranch, rightBranch)
    cosTheta = dotProduct / (mag(leftBranch) * mag(rightBranch))
    angle = math.acos(cosTheta)
    reflexAngle = (2 * math.pi) - angle
    halfReflex = (0.5) * reflexAngle
    tempPoint = rotate(leftBranch, halfReflex, -perpVector)
    normTempVector = tempPoint.norm()
    tempVecToBeAdded = (dist / 2) * normTempVector
    newVectors = []
    (vec1X, vec1Y, vec1Z) = (currAtomVec.x + perpVecToBeAdded.x + tempVecToBeAdded.x,
    currAtomVec.y + perpVecToBeAdded.y + tempVecToBeAdded.y, currAtomVec.z + perpVecToBeAdded.z + tempVecToBeAdded.z)
    (vec2X, vec2Y, vec2Z) = (currAtomVec.x - perpVecToBeAdded.x + tempVecToBeAdded.x, 
    currAtomVec.y - perpVecToBeAdded.y + tempVecToBeAdded.y, currAtomVec.z - perpVecToBeAdded.z + tempVecToBeAdded.z)
    newVectors.append(vector(vec1X, vec1Y, vec1Z)), newVectors.append(vector(vec2X, vec2Y, vec2Z))
    for i in range(len(atomToBeAdded)):
        color = figureColor(atomToBeAdded[i])
        newAtom = sphere(pos = newVectors[i], color = color)
        newAtom.label = atomToBeAdded[i]
    drawBonds(currAtomVec, newVectors, bondLst)


def popupMsg(msg):
    popup = Tk()
    popup.wm_title('Error!')
    label = Label(popup, text = msg, font = 'Arial 12 bold')
    label.pack(side = 'top', fill = 'x', pady = 10)
    B1 = Button(popup, text = 'Okay', command = popup.destroy)
    B1.pack()
    popup.mainloop()
    
#This is the function that controls all drawing
#Matrix is a dictionary of the following form:
# {Atom1 : [(Atom2, 'sb'), (Atom3, 'sb'), (Atom4, 'sb'), trigPlanar]...}
def mainFlow(matrix, filename, weight, dist = 10): #No rings
    clearScreen(filename, weight)
    placedAtoms = set() #Set of labels of placed atoms in the scene
    finishedAtoms = set() #Set of atoms whose bonding is entirely done
    counter, limit = 0, 1000
    while (set(matrix.keys()) != placedAtoms) and counter < limit: #Blocks until the entire molecule has been plotted
        for atom in matrix:
            centralAtom = atom
            connectingAtoms = matrix[atom]
            geometry = matrix[atom][-1] #Geometry is always the last element in the list
            atomLst, bondLst = unpackAtomBondList(connectingAtoms[:-1]) #atom list w/o central atom
            fullAtomLst = [centralAtom] + atomLst #list including the central atom
            if len(scene.objects) == 0: #The first atom on the scene
                if geometry == 'linear':
                    startingLinear(dist, fullAtomLst, bondLst)
                elif geometry == 'trigPlanar':
                    startingTrigPlanar(dist, fullAtomLst, bondLst)
                elif geometry == 'tetrahedral':
                    startingTetrahedral(dist, fullAtomLst, bondLst)
                elif geometry == 'trigBiPyramid':
                    startingTrigBiPyramid(dist, fullAtomLst, bondLst)
                elif geometry == 'octahedral':
                    startingOctahedral(dist, fullAtomLst, bondLst)
                elif geometry == None:
                    placingNone(fullAtomLst, dist, bondLst)
                finishedAtoms.add(centralAtom)
                for atom in fullAtomLst: placedAtoms.add(atom)
            else: #The case where we're building on
                if (atom in placedAtoms) and (atom not in finishedAtoms) and atomLst != []:
                    currAtom, prevAtom, badIndex = atom, None, None
                    for i in range(len(atomLst)):
                        label = atomLst[i]
                        if label in finishedAtoms:
                            prevAtom, badIndex = label, i
                    atomLst.pop(badIndex), bondLst.pop(badIndex)
                    #Add method here for finding existing bonds
                    if atomLst == []:
                        pass
                    elif geometry == 'linear':
                        addingLinear(currAtom, prevAtom, dist, atomLst, bondLst)
                    elif geometry == 'trigPlanar':
                        addingTrigPlanar(currAtom, prevAtom, dist, atomLst, bondLst)
                    elif geometry == 'tetrahedral':
                        addingTetrahedral(currAtom, prevAtom, dist, atomLst, bondLst)
                    elif geometry == 'trigBiPyramid':
                        addingTrigBiPyramid(currAtom, prevAtom, dist, atomLst, bondLst)
                    elif geometry == 'octahedral':
                        addingOctahedral(currAtom, prevAtom, dist, atomLst, bondLst)
                    finishedAtoms.add(currAtom)
                    for atom in atomLst: placedAtoms.add(atom)
                else: continue
        counter += 1
    if counter == limit:
        clearScreen('', 0)
        popupMsg('Your Molecule has errors!')

#Checks for existence in the scene
def checkExistence(ID):
    for obj in scene.objects:
        if isinstance(obj, sphere) and obj.label == ID: return True
    return False

def checkRingMembership(atomID, matrix):
    for ring in matrix['rings']:
        for ID in ring:
            if ID == atomID: return True
    return False

def checkIfFinished(matrix, atom):
    objLabelSet = set([obj.label for obj in scene.objects if isinstance(obj, sphere)])
    for i in range(len(matrix[atom])):
        value = matrix[atom][i]
        if type(value) == tuple:
            if value[0] not in objLabelSet: return False
    return True

#Code below here deals with rings
#Assumption: Rings constituent atoms have only tetrahedral geometry
def mainFlowWithRings(matrix, filename, weight, dist = 10):
    clearScreen(filename, weight)
    placedAtoms = set()
    finishedAtoms = set()
    counter, limit = 0, 1000
    keySet = set(matrix.keys())
    keySet.remove('rings')
    ringStarterIndex = 1
    geomSet = set()
    geomSet.add('linear'), geomSet.add('trigPlanar'), geomSet.add('tetrahedral'), geomSet.add('trigBiPyramid'), geomSet.add('octahedral')
    while (keySet != placedAtoms) and counter < limit:
        if len(scene.objects) == 0: #No previous atoms on the screen, we draw the first ring
            ring = matrix['rings'][0]
            startingRings(dist / 2, ring, matrix)
            for atoms in ring:
                placedAtoms.add(atoms)
                if checkIfFinished(matrix, atoms): #There are no atoms attached to it besides what's in the ring
                    finishedAtoms.add(atoms)
        else: 
            for ring in matrix['rings'][ringStarterIndex:]: #Check to see if there is a way to build a ring
                for atm in ring:
                    if checkExistence(atm) and atm not in finishedAtoms:
                        connectedAtomsToAtom = matrix[atm]
                        currAtom, prevAtom = atm, None
                        for i in range(len(connectedAtomsToAtom)): 
                            value = connectedAtomsToAtom[i]
                            if type(value) == tuple:
                                (atm2, bond) = value    
                                if atm2 in placedAtoms: prevAtom = atm2
                            elif type(value) == str: pass
                        ringLstCpy = copy(ring)
                        addingRings(currAtom, prevAtom, ringLstCpy, dist / 2) 
                        for atm in ringLstCpy:
                            placedAtoms.add(atm)
                            if checkIfFinished(matrix, atm):
                                finishedAtoms.add(atm)
                        ringStarterIndex += 1
                        break
            
            for obj in scene.objects: #Building atoms off of rings
                if isinstance(obj, sphere):
                    atmIdentif = obj.label
                    if (atmIdentif not in finishedAtoms) and type(matrix[atmIdentif][-1]) != str and checkRingMembership(atmIdentif, matrix):
                        atomsList = matrix[atmIdentif]
                        stringIndex = None
                        for i in range(len(atomsList)):
                            if type(atomsList[i]) == str:
                                stringIndex = i
                        ringID = atomsList[stringIndex] #E.g. 'ring0', 'ring1', etc.
                        ringIDNum = int(re.split('(\d+)', ringID)[1])
                        atomsToBeAttached = atomsList[stringIndex + 1 : ]
                        newAtomList, newBondList = unpackAtomBondList(atomsToBeAttached)
                        ringOfAtom = matrix['rings'][ringIDNum]
                        atmIDInRing = ringOfAtom.index(atmIdentif)
                        adjacentAtoms = [] #The one on the left and the one on the right will be in the adjacentAtoms
                        adjacentAtoms.append(ringOfAtom[atmIDInRing - 1]), adjacentAtoms.append(ringOfAtom[(atmIDInRing + 1) % len(ringOfAtom)])
                        if len(atomsList) - 1 - stringIndex == 1:
                            addingAtomsToRingsOne(adjacentAtoms, atmIdentif, newAtomList, dist, newBondList) 
                        elif len(atomsList) - 1 - stringIndex == 2:
                            addingAtomsToRingsTwo(adjacentAtoms, atmIdentif, newAtomList, dist, newBondList) 
                        finishedAtoms.add(atmIdentif)
                        for atoms in newAtomList:
                            placedAtoms.add(atoms)

                for atom in matrix: #Building atoms off atoms
                    if (type(matrix[atom][-1]) == str) and (matrix[atom][-1] in geomSet) and ((atom in placedAtoms) and (atom not in finishedAtoms)):
                        currAtom, prevAtom, badIndex = atom, None, None
                        atomListCurrent = matrix[atom]
                        geometry = matrix[atom][-1]
                        for i in range(len(atomListCurrent)):
                            value = atomListCurrent[i]
                            if type(value) == tuple:
                                (label, bond) = value
                                if label in finishedAtoms:
                                    prevAtom, badIndex = label, i
                        tempAtomList = copy(atomListCurrent)
                        tempAtomList.pop(badIndex), tempAtomList.pop()
                        unpackedAtomList, unpackedBondList = unpackAtomBondList(tempAtomList)
                        if len(tempAtomList) == 1 and type(tempAtomList[0]) == str:
                            pass
                        elif geometry == 'linear':
                            addingLinear(currAtom, prevAtom, dist, unpackedAtomList, unpackedBondList)
                        elif geometry == 'trigPlanar':
                            addingTrigPlanar(currAtom, prevAtom, dist, unpackedAtomList, unpackedBondList)
                        elif geometry == 'tetrahedral':
                            addingTetrahedral(currAtom, prevAtom, dist, unpackedAtomList, unpackedBondList)
                        elif geometry == 'trigBiPyramid':
                            addingTrigBiPyramid(currAtom, prevAtom, dist, unpackedAtomList, unpackedBondList)
                        elif geometry == 'octahedral':
                            addingOctahedral(currAtom, prevAtom, dist, unpackedAtomList, unpackedBondList)
                        finishedAtoms.add(currAtom)
                        for atom in unpackedAtomList: placedAtoms.add(atom)
        counter += 1
    if counter == limit:
        clearScreen('', 0)
        popupMsg('Your Moleule has errors!')


def convert2DNumpyArrayToList(numpyArray):
    result = list(numpyArray)
    for i in range(len(result)):
        result[i] = list(result[i])
    return result

#For finished Data
#CartesianMatrix of form [[Atom, X, Y, Z], [Atom2, X, Y, Z], ...]
#BondMatrix of form [[Atom1BT, BL], ...]
def mainFLowWithFinishedData(cartesianMatrix, bondMatrix, filename, weight):
    #plot all the atoms 
    clearScreen(filename, weight)
    cartesianMatrix = convert2DNumpyArrayToList(cartesianMatrix)
    bondMatrix = convert2DNumpyArrayToList(bondMatrix)
    for atomData in cartesianMatrix:
        atomLabel = atomData[0]
        position = vector(atomData[1], atomData[2], atomData[3])
        color = figureColor(atomLabel)
        newAtom = sphere(pos = position, color = color)
        newAtom.label = atomLabel
    for bond in bondMatrix:
        bondMetaDat, bondLen = bond[0], bond[1]
        bondMetaDatUnpacked = re.split('(\d+)', bondMetaDat)
        bondType = bondMetaDatUnpacked.pop()
        atomList = []
        for i in range(len(bondMetaDatUnpacked)):
            if i % 2 == 0:
                atomAndNum = bondMetaDatUnpacked[i] + bondMetaDatUnpacked[i + 1]
                atomList.append(atomAndNum)
        centAtom, connectedAtom = atomList[0], atomList[1]
        centAtomPos, connectedAtomPos = None, None
        for obj in scene.objects:
            if isinstance(obj, sphere):
                if obj.label == centAtom:
                    centAtomPos = obj.pos
                elif obj.label == connectedAtom:
                    connectedAtomPos = obj.pos
        drawBonds(centAtomPos, [connectedAtomPos], [bondType])



    
    
    