# -*- coding: utf-8 -*-
#File that handles all the lewis structure drawings and most of the user interface

import Geoms
from tkinter import *
import math
from mendeleev import element
import copy
import re 
import pickle
import os
import os.path
from tkinter import filedialog
from targetOutputBaseNew import runGetAnglesAndBondsBase
import numpy as np
from PIL import Image, ImageTk


####################################
#Classes and Methods
####################################

#Classes for atom and bonds + subclasses of bonds

class Atom:
    #Model
    #An atom has the atomic symbol, its number, its center, and its radius + text scale factor
    def __init__(self, atom, cX, cY, radius, ID, color = 'black'):
        #Drawing data
        self.atom = atom
        self.cX, self.cY = cX, cY
        self.r = radius
        self.color = color
        #Valency and geometry data
        self.ID = ID #ID of the form element symbol + number
        self.numSingleBonds = 0
        self.numDoubleBonds = 0
        self.numTripleBonds = 0
        self.lonePairs = []
        self.cardinals = Atom.getCardinals(self.cX, self.cY, self.r)
        self.valencyLimit, self.bondLimit = Atom.getValencyData(self.atom)
        self.takenPoints = set()
        self.bondDomainNumber = (self.numDoubleBonds + self.numSingleBonds + self.numTripleBonds)
        self.currValency = 0
        self.currBonds = 0
        self.lpSingles = self.bondLimit
        self.lpDoubles = self.valencyLimit - self.bondLimit
        if self.valencyLimit == self.bondLimit == 6:
            self.lpSingles = 0
            self.lpDoubles = self.valencyLimit
        #Steric number and geometry data
        self.stericNumber = self.bondDomainNumber + len(self.lonePairs) #Num of bonded atoms + num lone pairs on this atom
        self.geometry = Atom.getGeometryFromSterics(self.stericNumber)
    
    #Generates a list of contact points for the atom, just once at the beginning
    #   Each atom has 8 contact points to start

    #View
    def draw(self, canvas):
        canvas.create_oval(self.cX - self.r, self.cY - self.r, self.cX + self.r, 
                           self.cY + self.r, width = 0, fill = 'white')
        canvas.create_text(self.cX, self.cY, font = 'Arial ' + 
                           str(int(self.r)) + ' bold', fill = self.color, text = self.atom)
        Atom.drawLonePairs(self, canvas)
    
    #Controller

    def generateLonePair(self, data):
        badPoints = [] #A tuple of the bad points that links a bond to the atom
        coordsToCompare = None
        for bond in data.bondList:
            if self in bond.connectingAtoms:
                #If we dealing with double or triple bond, we just want the coords of the center line
                if isinstance(bond, DoubleBond) or isinstance(bond, TripleBond): 
                    coordsToCompare = bond.coordTupleSet3
                #If we are dealing with only a single bond
                elif isinstance(bond, Bond) and not (isinstance(bond, DoubleBond) or isinstance(bond, TripleBond)):
                    coordsToCompare = bond.coordTupleSet
                if isinstance(bond, DoubleBond) or isinstance(bond, TripleBond):
                    for (x, y) in coordsToCompare: 
                        if almostEqual(distance(self.cX, self.cY, x, y), self.r):
                            badPoints.append((x, y))
                elif isinstance(bond, Bond) and not (isinstance(bond, DoubleBond) or isinstance(bond, TripleBond)):
                        badPoints.append(Atom.getEdgeLocationsSingleBond(coordsToCompare, self)) 
        #Drawing on the cardinals                       
        if len(badPoints) == 0 and self.cardinals != []:
            for (x, y) in self.cardinals:
                if (x, y) not in self.takenPoints:
                    if self.currBonds + len(self.lonePairs) < self.lpSingles:
                        self.lonePairs.append((x, y, 's'))
                    elif self.currBonds == self.lpSingles or (len(self.lonePairs)) == self.lpSingles or self.currBonds + len(self.lonePairs) == self.lpSingles:
                        self.lonePairs.append((x, y, 'd'))
                    self.takenPoints.add((x, y))                  
                    self.cardinals = self.cardinals[1:]
                    return
        elif (len(badPoints) == 0 and self.cardinals == []) or len(badPoints) != 0:
            winningPoint, maximumDistance = None, 0
            tempDistance = 0
            for sector in range(12):
                angle = sector * math.pi / 6
                newX, newY = self.cX + self.r * math.cos(angle), self.cY - self.r * math.sin(angle)
                if (newX, newY) not in self.takenPoints:
                    if len(badPoints) != 0:
                        for (x, y) in badPoints:
                            tempDistance += distance(newX, newY, x, y)
                        if tempDistance > maximumDistance:
                            maximumDistance, winningPoint = tempDistance, (newX, newY)   
                    elif len(badPoints) == 0:
                        winningPoint = (newX, newY)
                tempDistance = 0
            if self.currBonds + len(self.lonePairs) < self.lpSingles:
                self.lonePairs.append((winningPoint[0], winningPoint[1], 's'))
            elif self.currBonds >= self.lpSingles or (len(self.lonePairs)) >= self.lpSingles or self.currBonds + len(self.lonePairs) >= self.lpSingles:
                self.lonePairs.append((winningPoint[0], winningPoint[1], 'd'))
            self.takenPoints.add(winningPoint)
            
            
            
    @staticmethod
    def getCardinals(x, y, r):
        north, south = (x, y - r), (x, y + r)
        west, east = (x - r, y), (x + r, y)
        return [north, south, west, east]
    
    @staticmethod
    def drawLonePairs(atom, canvas):
        buffer = 2
        for (xCoord, yCoord, lpType) in atom.lonePairs:
            if lpType == 'd':
                vector = [xCoord - atom.cX, yCoord - atom.cY]
                normalVector = vectorNormalize(vector)
                leftNorVec, rightNorVec = computerOrthogonals(normalVector)
                firstPointX, firstPointY = xCoord + buffer * leftNorVec[0], yCoord + buffer * leftNorVec[1]
                secondPointX, secondPointY = xCoord + buffer * rightNorVec[0], yCoord + buffer * rightNorVec[1]
                canvas.create_oval(firstPointX - 1, firstPointY - 1, firstPointX + 1, firstPointY + 1, fill = 'black')
                canvas.create_oval(secondPointX - 1, secondPointY - 1, secondPointX + 1, secondPointY + 1, fill = 'black')
            elif lpType == 's':
                canvas.create_oval(xCoord - 1, yCoord - 1, xCoord + 1, yCoord + 1, fill = 'black')
        
    @staticmethod
    def getEdgeLocationsSingleBond(coordSet, atom):
        startingPoint = None
        for (x, y) in coordSet:
            if math.isclose(atom.cX, x) and math.isclose(atom.cY, y):
                startingPoint = (x, y)
        endPoint = None
        for (x, y) in coordSet:
            if (x, y) != startingPoint: endPoint = (x, y)
        vector = [endPoint[0] - startingPoint[0], endPoint[1] - startingPoint[1]]
        normalVector = vectorNormalize(vector)
        badX = atom.cX + atom.r * normalVector[0]
        badY = atom.cY + atom.r * normalVector[1]
        return (badX, badY)
        
    @staticmethod
    def getValencyData(elem):
        elemObj = element(elem)
        valencyLimit, bondLimit = None, None
        #Check what period the element is in. Greater than 3, expanded valency up to 
        #   octahedral
        if elemObj.period >= 3:
            valencyLimit, bondLimit = 6, 6
        #Elements in period 2 require 
        elif elemObj.period == 2:
            groupName = elemObj.group.name
            if groupName == 'Boron group': valencyLimit, bondLimit = 3, 3
            elif groupName == 'Carbon group': valencyLimit, bondLimit = 4, 4
            elif groupName == 'Pnictogens': valencyLimit, bondLimit = 4, 3
            elif groupName == 'Chalcogens': valencyLimit, bondLimit = 4, 2
            elif groupName == 'Halogens': valencyLimit, bondLimit = 4, 1
        #If less than two, the only one we care about is hydrogen
        elif elemObj.period == 1:
            valencyLimit, bondLimit = 1, 1
        return valencyLimit, bondLimit
    
    @staticmethod
    #Gets the molecular geometry from the steric number and bondDomain number (number of bonds)
    def getGeometryFromSterics(steric):
        if steric == 2: return 'linear'
        elif steric == 3: return 'trigPlanar'
        elif steric == 4: return 'tetrahedral'
        elif steric == 5: return 'trigBiPyramid'
        elif steric == 6: return 'octahedral'
    
#The parent bond class is the same as a single bond
class Bond:
    #Model
    def __init__(self, coordSet, atom1, atom2):
        self.x1, self.y1, self.x2, self.y2 = coordSet
        self.coordTupleSet = set()
        self.coordTupleSet.add((self.x1, self.y1))
        self.coordTupleSet.add((self.x2, self.y2))
        #Keep track of the atoms that the bond is linking
        self.connectingAtoms = [atom1, atom2]
    
    #View
    #Will update this method later with how to draw double and triple bonds
    def draw(self, canvas):
        canvas.create_line(self.x1, self.y1, self.x2, self.y2, width = 2)

class DoubleBond(Bond):
    def __init__(self, coordSet, coordSet2, coordSet3, atom1, atom2):
        super().__init__(coordSet, atom1, atom2)
        self.secX1, self.secY1, self.secX2, self.secY2 = coordSet2
        self.coordTupleSet2 = set()
        self.coordTupleSet2.add((self.secX1, self.secY1))
        self.coordTupleSet2.add((self.secX2, self.secY2))
        self.thrX1, self.thrY1, self.thrX2, self.thrY2 = coordSet3
        self.coordTupleSet3 = set()
        self.coordTupleSet3.add((self.thrX1, self.thrY1))
        self.coordTupleSet3.add((self.thrX2, self.thrY2))
        pass
    
    def draw(self, canvas):
        canvas.create_line(self.x1, self.y1, self.x2, self.y2, width = 2)
        canvas.create_line(self.secX1, self.secY1, self.secX2, self.secY2, 
                           width = 2)
    
class TripleBond(DoubleBond):
    def __init__(self, coordSet, coordSet2, coordSet3, atom1, atom2):
        super().__init__(coordSet, coordSet2, coordSet3, atom1, atom2)
        pass
    
    def draw(self, canvas):
        super().draw(canvas)
        canvas.create_line(self.thrX1, self.thrY1, self.thrX2, self.thrY2, 
                           width = 2)

class Ring:
    def __init__(self, cX, cY, vertices, radius = 80):
        self.cX, self.cY, self.vertices, self.r = cX, cY, vertices, radius
        self.vertexList = []
        self.vertexList = Ring.generatePolygonCoordinates(self.cX, self.cY, self.vertices, self.r)
        self.vertexIDList = [i for i in range(len(self.vertexList))]
        self.connectedAtoms = [None for i in range(len(self.vertexList))] #A series of nonetypes that we can replace gradually
    
    def draw(self, canvas):
        vertexListUnpacked = []
        for (x, y) in self.vertexList:
            vertexListUnpacked.append(x), vertexListUnpacked.append(y)
        canvas.create_polygon(vertexListUnpacked, fill = '', outline = 'black', width = 2)
            
    @staticmethod
    def generatePolygonCoordinates(x, y, vertices, radius):
        vertexList = []
        sectorAngle = (2 * math.pi) / vertices
        for i in range(vertices):
            angle = i * sectorAngle
            newX, newY = x + (radius * math.cos(angle)), y - (radius * math.sin(angle))
            vertexList.append((newX, newY))
        return vertexList

# =============================================================================
# Main program functions
# =============================================================================
    
#Comparing if two values are around equal based on 5% of max value
def almostEqual(x, y):
    maxValue = max(x, y)
    return abs(x - y) < 0.08 * maxValue

#Handles comparisons near zero, 
def almostEqualNearZero(x, y, data):
    return abs(x - y) < 0.10 * (2 * data.atomRadius)
        
#Euclidian distance
def distance(x1, y1, x2, y2):
    return (((x2 - x1) ** 2) + ((y2 - y1) ** 2)) ** 0.5

#m = (y2 - y1) / (x2 - x1)
def calculateSlope(x1, y1, x2, y2, data):
    if almostEqualNearZero((x2 - x1), 0, data):
        return 'Vertical'
    else: return (y2 - y1) / (x2 - x1)

#y = mx + b, finds b
def findIntercept(slope, x, y):
    return y - (slope * x)


#Function for normalizing vectors
def vectorNormalize(vector):
    vecMagnitude = (((vector[0]) ** 2) + ((vector[1]) ** 2)) ** 0.5
    normalVector = copy.copy(vector)
    for i in range(len(normalVector)):
        normalVector[i] /= vecMagnitude
    return normalVector

#Function for giving anti-parallel vectors
def vectorReverse(vector):
    reverseVector = copy.deepcopy(vector)
    for i in range(len(reverseVector)):
        reverseVector[i] *= -1
    return reverseVector

#In 2D space, finds the two vectors that are perpendicular to a parent vector
def computerOrthogonals(vector):
    #There are two possible perpendicular vectors, and we will find both by 
    #   switching and negating the original vector's x + y components. The original
    #   vector is <x, y>. The left will be <y, -x> and the right will be <-y, x>
    leftPerpVector, rightPerpVector = copy.copy(vector), copy.copy(vector)
    #computing the left vector
    leftPerpVector[0], leftPerpVector[1] = leftPerpVector[1], leftPerpVector[0]
    leftPerpVector[1] = -leftPerpVector[1]
    rightPerpVector[0], rightPerpVector[1] = rightPerpVector[1], rightPerpVector[0]
    rightPerpVector[0] = -rightPerpVector[0]
    return leftPerpVector, rightPerpVector

#Function that generates the coordinate sets for bonds
def generateCoordSets(data, flag):
    if flag == 'db':
        sepFac = 2
    elif flag == 'tb':
        sepFac = 5
    atom1, atom2 = data.atoms[data.selectedAtoms[0]], data.atoms[data.selectedAtoms[1]]
    cX1, cY1, cX2, cY2 = atom1.cX, atom1.cY, atom2.cX, atom2.cY
    #Vector between the two centers of two atoms, pointed in the direction
    #   of the second atom
    vector = [cX2 - cX1, cY2 - cY1]
    #Returns normal vector of magnitude 1
    normalVector = vectorNormalize(vector)
    #Anti-parallel of normalVector, pointing towards first atom
    antiNormalVector = vectorReverse(normalVector)
    stCenterX, stCenterY = cX1 + data.atomRadius * normalVector[0],\
    cY1 + data.atomRadius * normalVector[1]
    endCenterX, endCenterY = cX2 + data.atomRadius * antiNormalVector[0],\
    cY2 + data.atomRadius * antiNormalVector[1]
    #The imaginary line b/w the two bonding lines is thought of as 
    #   (stCenterX, stCenterY) --> (endCenterX, endCenterY)
    leftVec, rightVec = computerOrthogonals(normalVector)
    firstX1, firstY1 = stCenterX + sepFac * leftVec[0], stCenterY + sepFac * leftVec[1]
    firstX2, firstY2 = endCenterX + sepFac * leftVec[0], endCenterY + sepFac * leftVec[1]
    secondX1, secondY1 = stCenterX + sepFac * rightVec[0], stCenterY + sepFac * rightVec[1]
    secondX2, secondY2 = endCenterX + sepFac * rightVec[0], endCenterY + sepFac * rightVec[1]
    coordSet = (firstX1, firstY1, firstX2, firstY2)
    coordSet2 = (secondX1, secondY1, secondX2, secondY2)
    #CoordSet for the center line is going to be the third coordSet
    coordSet3 = (stCenterX, stCenterY, endCenterX, endCenterY)
    return coordSet, coordSet2, coordSet3, atom1, atom2


#If you remove an atom, we also want to remove all the bonds tied to that 
#   atom, and decrease the counts of all other connected atoms accordingly
def removeBadBondsAtomDependent(data, target):
    remAtom = data.atoms[target]
    badBonds = set()
    for i in range(len(data.bondList)):
        bond = data.bondList[i]
        if remAtom in bond.connectingAtoms:
            badBonds.add(i)
    newBondList = []
    for i in range(len(data.bondList)):
        bond = data.bondList[i]
        if i not in badBonds: newBondList.append(bond)
        else:
            for atom in bond.connectingAtoms:
                if isinstance(bond, Bond) and not  (isinstance(bond, DoubleBond) or isinstance(bond, TripleBond)):
                    atom.numSingleBonds -= 1
                elif isinstance(bond, DoubleBond) and not isinstance(bond, TripleBond):
                    atom.numDoubleBonds -= 1
                elif isinstance(bond, TripleBond):
                    atom.numTripleBonds -= 1
                updateCurrValency(atom)
    data.bondList = newBondList

def singleBondCheck(mousex, mousey, bond, data):
    farthestRight = max(bond.x1, bond.x2)
    farthestLeft = min(bond.x1, bond.x2)
    if mousex > farthestRight or mousex < farthestLeft: pass
    else:
        slope = calculateSlope(bond.x1, bond.y1, bond.x2, bond.y2, data)
        #Dealing with near vertical bonds
        if slope == 'Vertical':
            farthestUp = min(bond.y1, bond.y2)
            farthestDown = max(bond.y1, bond.y2)
            if farthestUp < mousey < farthestDown:
                return True
        else:
            intCept = findIntercept(slope, bond.x1, bond.y1)
            eqOutput = (slope * mousex) + intCept
            #If the mouse position is close enough to the bond position
            #   for math.isclose, the rel_tol is a percentage
            if almostEqual(eqOutput, mousey):
                return True
            
def doubleTripleCheck(mousex, mousey, bond, data):
    firstX1, firstY1, firstX2, firstY2 = bond.x1, bond.y1, bond.x2, bond.y2
    secX1, secY1, secX2, secY2 = bond.secX1, bond.secY1, bond.secX2, bond.secY2
    thrX1, thrY1, thrX2, thrY2 = bond.thrX1, bond.thrY1, bond.thrX2, bond.thrY2
    farthestRight = max(firstX1, firstX2, secX1, secX2, thrX1, thrX2)
    farthestLeft = min(firstX1, firstX2, secX1, secX2, thrX1, thrX2)
    if mousex > farthestRight or mousex < farthestLeft: pass
    else:   
        slope = calculateSlope(firstX1, firstY1, firstX2, firstY2, data)
        if slope == 'Vertical':
            farthestUp = min(firstY1, firstY2, secY1, secY2, thrY1, thrY2)
            farthestDown = max(firstY1, firstY2, secY1, secY2, thrY1, thrY2)
            if farthestUp < mousey < farthestDown:
                return True
        else:
            #Wildly unequal intercepts
            intCept1 = findIntercept(slope, firstX1, firstY1)
            intCept2 = findIntercept(slope, secX1, secY1)
            intCept3 = findIntercept(slope, thrX1, thrY1)
            eqOutput1 = (slope * mousex) + intCept1
            eqOutput2 = (slope * mousex) + intCept2
            eqOutput3 = (slope * mousex) + intCept3
            maxOut = max(eqOutput1, eqOutput2, eqOutput3)
            minOut = min(eqOutput1, eqOutput2, eqOutput3)
            proxCondition = (almostEqual(eqOutput1, mousey) or almostEqual(eqOutput2, mousey)\
            or almostEqual(eqOutput3, mousey)) or (minOut < mousey < maxOut)
            if proxCondition:
                return True

#Removing bonds just by clicking on them, checks slope formulas. Takes in event.x,
#   event.y, and data
def removeBadBondsNonAtomDependent(mousex, mousey, data, flag):
    badBondIndSet = set()
    for i in range(len(data.bondList)):
        bond = data.bondList[i]
        #if we are dealing with a single bond
        if isinstance(bond, Bond) and not (isinstance(bond, DoubleBond) or isinstance(bond, TripleBond)):
            if singleBondCheck(mousex, mousey, bond, data):
                if len(badBondIndSet) < 1:
                    badBondIndSet.add(i)
                    flag.append(True)
        #Double bonds and triple bonds require that we check three lines
        elif isinstance(bond, DoubleBond) or isinstance(bond, TripleBond):
            if doubleTripleCheck(mousex, mousey, bond, data):
                if len(badBondIndSet) < 1:
                    badBondIndSet.add(i)
                    flag.append(True)
    newBondList = []
    for i in range(len(data.bondList)):
        bond = data.bondList[i]
        if i not in badBondIndSet: newBondList.append(bond)
        else:
            for atom in bond.connectingAtoms:
                if isinstance(bond, Bond) and not (isinstance(bond, DoubleBond) or isinstance(bond, TripleBond)):
                    atom.numSingleBonds -= 1
                elif isinstance(bond, DoubleBond) and not isinstance(bond, TripleBond):
                    atom.numDoubleBonds -= 1
                elif isinstance(bond, TripleBond):
                    atom.numTripleBonds -= 1
                updateCurrValency(atom)
    data.bondList = newBondList

#TODO: Finish this method
def removeRings(mouseX, mouseY, data):
    badRingIndex = None
    for i in range(len(data.rings)):
        ring = data.rings[i]
        cX, cY = ring.cX, ring.cY
        dist = distance(mouseX, mouseY, cX, cY)
        if dist < data.atomRadius:
            badRingIndex = i
            break
    if badRingIndex == None:
        return
    else:
        badRing = data.rings[badRingIndex]
        for atom in badRing.connectedAtoms:
            if atom != None:
                atom.numSingleBonds -= 2
                updateCurrValency(atom)
        data.rings.pop(badRingIndex)

#Initializes data for the program
def init(data):
    data.atoms = []
    data.atom = 'C' 
    data.cursorMode = 'atomPlacement'
    data.atomRadius = 20
    data.selectedAtoms = []
    data.bondType = 'single'
    data.bondList = []
    data.rightCursorMode = 'rmAtom'
    data.atomDict = dict()
    data.rings = []
    data.ringVertices = 3
    data.moleculeFormula = ""
    data.moleculeWeight = 0
    data.order = ['C', 'H', 'N', 'O']
    data.stringVar = StringVar()
    data.stringVar.set("Formula: %s    Weight: %0.2f g/mol    Cursor Mode: %s    Atom: %s" % 
                       (data.moleculeFormula, data.moleculeWeight, data.cursorMode, data.atom))

def partialInit(data):
    data.atoms = []
    data.atom = 'C' 
    data.cursorMode = 'atomPlacement'
    data.selectedAtoms = []
    data.bondType = 'single'
    data.bondList = []
    data.atomDict = dict()
    data.rings = []
    data.ringVertices = 3
    data.moleculeFormula = ""
    data.moleculeWeight = 0
    data.stringVar.set("Formula: %s    Weight: %0.2f g/mol    Cursor Mode: %s    Atom: %s" % 
                       (data.moleculeFormula, data.moleculeWeight, data.cursorMode, data.atom))
    
def findNearestRingVertex(x, y, data):
    ringID, closestVertexID, minDistance = None, None, None
    for i in range(len(data.rings)): #The ring ID in data.rings
        ring = data.rings[i]
        for j in range(len(ring.vertexList)): #The vertexID inside the ring's vertex list
            (x2, y2) = ring.vertexList[j]
            dist = distance(x, y, x2, y2)
            if minDistance == None:
                minDistance, closestVertexID, ringID = dist, j, i
            elif minDistance != None:
                if dist < minDistance:
                    minDistance, closestVertexID, ringID = dist, j, i
    return ringID, closestVertexID, minDistance
    

#Finds the closest atom to a mouseclick
def findClosestAtom(x, y, data):
    minDistance = None
    closestAtomIndex = None
    for i in range(len(data.atoms)):
        atom = data.atoms[i]
        x2, y2 = atom.cX, atom.cY
        dist = distance(x, y, x2, y2)
        if minDistance == None:
            minDistance, closestAtomIndex = dist, i
        elif minDistance != None:
            if dist < minDistance: 
                minDistance, closestAtomIndex = dist, i
    return minDistance, closestAtomIndex


#Checks to see if the atom to be placed overlaps with any other atom
def checkOverlap(atom, data):
    #Simple check to see if the atoms overlap
    x1, y1 = atom.cX, atom.cY
    for atom in data.atoms:
        x2, y2 = atom.cX, atom.cY
        #Overlap is found
        if distance(x1, y1, x2, y2) <= 2 * data.atomRadius:
            return True
    return False

#Checks to see if the atom to be placed overlaps any existing bonds, based off
#   a series of points around the atom
def checkBondOverlap(atom, data):
    cX, cY = atom.cX, atom.cY
    for bond in data.bondList:
        if isinstance(bond, Bond) and not ((isinstance(bond, DoubleBond) or isinstance(bond, TripleBond))):
            #Chose 12 b/c 2pi / (pi/6) = 12
            for sector in range(12):
                angle = sector * math.pi / 6
                x, y = cX + atom.r * math.cos(angle), cY - atom.r * math.sin(angle)
                if singleBondCheck(x, y, bond, data):
                    return True
        elif isinstance(bond, DoubleBond) or isinstance(bond, TripleBond):
            for sector in range(12):
                angle = sector * math.pi / 6
                x, y = cX + atom.r * math.cos(angle), cY - atom.r * math.sin(angle)
                if doubleTripleCheck(x, y, bond, data):
                    return True
    return False
                    
def updateBondNumber(data): #Updates the number of bonds that an atom has attached to it
    selectedAtomIDs = [data.atoms[data.selectedAtoms[0]].ID, data.atoms[data.selectedAtoms[1]].ID]
    selectedAtomIDSet = set(selectedAtomIDs)
    for bond in data.bondList:
        joinedAtomIDs = [atom.ID for atom in bond.connectingAtoms]
        joinedAtomIDSet = set(joinedAtomIDs)
        if joinedAtomIDSet == selectedAtomIDSet:
            for atom in data.atoms:
                if atom.ID in selectedAtomIDSet:
                    if isinstance(bond, Bond) and not ((isinstance(bond, DoubleBond) or isinstance(bond, TripleBond))):    
                        atom.numSingleBonds += 1
                    elif isinstance(bond, DoubleBond) and not isinstance(bond, TripleBond):
                        atom.numDoubleBonds += 1
                    elif isinstance(bond, TripleBond):
                        atom.numTripleBonds += 1
                    updateCurrValency(atom)
                        
#Updates the current valency and steric number and others
def updateCurrValency(atom):
    atom.currValency = (atom.numSingleBonds) + (2 * atom.numDoubleBonds) + (3 * atom.numTripleBonds) + len(atom.lonePairs)
    atom.currBonds = (atom.numSingleBonds) + (2 * atom.numDoubleBonds) + (3 * atom.numTripleBonds)
    atom.bondDomainNumber = atom.numSingleBonds + atom.numDoubleBonds + atom.numTripleBonds
    atom.stericNumber = atom.bondDomainNumber + len(atom.lonePairs)
    atom.geometry = Atom.getGeometryFromSterics(atom.stericNumber)

#TODO: Continue from here
def updateFormulaAndWeight(data):
    countDict = dict()
    formula, weight = '', 0
    #Gives us a count of all the elements
    for atom in data.atoms:
        elem = atom.atom
        if elem in countDict:
            countDict[elem] += 1
        elif elem not in countDict:
            countDict[elem] = 1
    keySet = set(countDict.keys())
    finishedAtoms = set()
    #Adds all the elements from data.order first
    for atomSym in data.order:
        if atomSym in keySet:    
            tempAddition = atomSym + str(countDict[atomSym])
            formula += tempAddition
            elemObj, count = element(atomSym), countDict[atomSym]
            weight += (count * elemObj.atomic_weight)
            finishedAtoms.add(atomSym)
    for atomos in countDict:
        if atomos not in finishedAtoms:
            tempAddition = atomos + str(countDict[atomos])
            formula += tempAddition
            elemObj, count = element(atomos), countDict[atomos]
            weight += (count * elemObj.atomic_weight)
            finishedAtoms.add(atomos)
    data.moleculeFormula, data.moleculeWeight = formula, weight
    data.stringVar.set("Formula: %s    Weight: %0.2f g/mol    Cursor Mode: %s    Atom: %s" % 
                       (data.moleculeFormula, data.moleculeWeight, data.cursorMode, data.atom))
        
#Generating popup messages
def popupMsg(msg):
    popup = Tk()
    #Icon taken from https://www.flaticon.com/free-icon/molecule_195742
    popup.iconbitmap(os.path.join(os.getcwd(), 'BackgroundImages', 'molecule.ico'))
    popup.resizable(width = False, height = False)
    popup.wm_title('Error!')
    label = Label(popup, text = msg, font = 'Arial 12 bold')
    label.pack(side = 'top', fill = 'x', pady = 10)
    B1 = Button(popup, text = 'Okay', command = popup.destroy)
    B1.pack()
    popup.mainloop()

#PAY CLOSE ATTENTION TO THIS FUNCTION
def mousePressed(event, data):
    # use event.x and event.y
    #Updates the number of the atom that is shown on screen
    if data.cursorMode == 'atomPlacement':
        element = data.atom 
        if element not in data.atomDict: data.atomDict[element] = 1
        else: data.atomDict[element] += 1
        ID = element + str(data.atomDict[element]) #Add in a unique identifier for each atom
        ringID, closestVertexID, minDistanceVertex = findNearestRingVertex(event.x, event.y, data)
        if minDistanceVertex != None and minDistanceVertex < data.atomRadius:
            atmX, atmY = data.rings[ringID].vertexList[closestVertexID]
            atom = Atom(element, atmX, atmY, data.atomRadius, ID)
            data.rings[ringID].connectedAtoms[closestVertexID] = atom
            data.atoms.append(atom)
            updateValencyInRing(data)
            updateFormulaAndWeight(data)
        else:
            atom = Atom(element, event.x, event.y, data.atomRadius, ID)
            if not checkOverlap(atom, data) and not checkBondOverlap(atom, data):
                data.atoms.append(atom)
                updateFormulaAndWeight(data)
    elif data.cursorMode == 'lonePairPlacement': #Placing lone pairs
        minDistance, closestAtomIndex = findClosestAtom(event.x, event.y, data)
        if minDistance == None: pass
        elif minDistance != None:
            if minDistance < data.atomRadius:
                atom = data.atoms[closestAtomIndex]
                if atom.currValency < atom.valencyLimit:
                    atom.generateLonePair(data)
                    updateCurrValency(atom)
                else: popupMsg('Maximum Valency Exceeded!')
    elif data.cursorMode == 'bondPlacement': #Creating bonds
        minDistance, closestAtomIndex = findClosestAtom(event.x, event.y, data)
        #if minDistance == None, no atom was found
        if minDistance == None: pass
        elif minDistance != None:
            if minDistance < data.atomRadius:
                atom = data.atoms[closestAtomIndex]
                atom.color = 'green'
                data.selectedAtoms.append(closestAtomIndex)
            else:
                for atom in data.atoms: atom.color = 'black' #clicking outside the molecule resets all atom colors to black
                data.selectedAtoms = []
                return 
            if len(data.selectedAtoms) == 2: #two atoms have been selected
                checkAndGenerateBonds(data)
            elif len(data.selectedAtoms) > 2: #more than two atoms have been selected
                data.selectedAtoms = [data.selectedAtoms[-1]] #Keep only the third one
                for i in range(len(data.atoms)):
                    if i not in data.selectedAtoms:
                        #Reset all colors to black
                        data.atoms[i].color = 'black'
    elif data.cursorMode == 'ringPlacement':
        ring = Ring(event.x, event.y, data.ringVertices)
        data.rings.append(ring)
        
def checkAndGenerateBonds(data):
    selectedAtomIDs = [data.atoms[data.selectedAtoms[0]].ID, data.atoms[data.selectedAtoms[1]].ID]
    selectedAtomIDSet = set(selectedAtomIDs)
    for bond in data.bondList:
        joinedAtomIDs = [atom.ID for atom in bond.connectingAtoms]
        joinedAtomIDSet = set(joinedAtomIDs)
        if selectedAtomIDSet == joinedAtomIDSet: return #A bond already exists between these two atoms
    for index in data.selectedAtoms:
        atom = data.atoms[index]
        if data.bondType == 'single':
            if not atom.currBonds <= (atom.bondLimit - 1):
                popupMsg('Maximum Valency Exceeded!')
                return
        elif data.bondType == 'double':
            if not atom.currBonds <= (atom.bondLimit - 2):
                popupMsg('Maximum Valency Exceeded!')
                return
        elif data.bondType == 'triple':
            if not atom.currBonds <= (atom.bondLimit - 3):
                popupMsg('Maximum Valency Exceeded!')
                return
    generateBond(data)
    updateBondNumber(data)

#General method for removing elements from the screen
def mouseRightClick(event, data):
    minDistance, closestAtomIndex = findClosestAtom(event.x, event.y, data)
    #if minDistance == None, no atom was found
    if minDistance == None: pass
    elif minDistance != None:
        if minDistance < data.atomRadius:
            removeBadBondsAtomDependent(data, closestAtomIndex) #FIX here so that all the atoms after the one that's popped are lowered in number by one
            if closestAtomIndex != len(data.atoms) - 1:
                closestAtomElement = data.atoms[closestAtomIndex].atom
                finalIndex = None
                for i in range(closestAtomIndex + 1, len(data.atoms)):
                    atom = data.atoms[i]
                    if atom.atom == closestAtomElement:
                        idStr = atom.ID
                        splitList = re.split('(\d+)', idStr)
                        number = int(splitList[1])
                        number -= 1
                        newID = splitList[0] + str(number)
                        atom.ID = newID
                        finalIndex = i
                #This was the last element of that kind in the list
                if finalIndex == None:
                    removedAtom = data.atoms[closestAtomIndex]
                    removeAtomFromRingData(removedAtom, data)
                    splitLst = re.split('(\d+)', removedAtom.ID)
                    number = int(splitLst[1])
                    data.atomDict[removedAtom.atom] = number - 1
                    data.atoms.pop(closestAtomIndex)
                #Another element of that kind was found later in the list
                elif finalIndex != None:
                    finalAtom = data.atoms[finalIndex]
                    #Handle the ring data first; the finalID is to deal with other atoms
                    removedAtom = data.atoms[closestAtomIndex]
                    removeAtomFromRingData(removedAtom, data)
                    finalID = finalAtom.ID
                    splitListFinal = re.split('(\d+)', finalID)
                    data.atomDict[finalAtom.atom] = int(splitListFinal[1])
                    data.atoms.pop(closestAtomIndex)
            else: #In the event that the atom you removed is the last one in the list
                removedAtom = data.atoms.pop(closestAtomIndex)
                removeAtomFromRingData(removedAtom, data)
                splitLst = re.split('(\d+)', removedAtom.ID)
                data.atomDict[removedAtom.atom] = int(splitLst[1]) - 1
            updateFormulaAndWeight(data)
            return
    #Here, we make it so that you can remove the bonds by right clicking on 
    flag = []
    removeBadBondsNonAtomDependent(event.x, event.y, data, flag)
    if flag == [True]:    
        updateFormulaAndWeight(data)
    else:
        removeRings(event.x, event.y, data)
        updateFormulaAndWeight(data)

        

#FIx this function
def removeAtomFromRingData(removedAtom, data):
    if len(data.rings) == 0:
        pass
    else:
        ringID, vertexID = None, None
        x, y = removedAtom.cX, removedAtom.cY
        for i in range(len(data.rings)):
            ring = data.rings[i]
            for j in range(len(ring.vertexList)):
                (x2, y2) = ring.vertexList[j]
                if x2 == x and y2 == y:
                    ringID, vertexID = i, j
        if ringID != None and vertexID != None:
            ring = data.rings[ringID]
            ring.connectedAtoms[vertexID] = None
            
def generateBond(data):
    #First, find the closest two contact points between atoms
    if data.bondType == 'single':
        atom1Index, atom2Index = data.selectedAtoms[0], data.selectedAtoms[1]
        atom1, atom2 = data.atoms[atom1Index], data.atoms[atom2Index]
        coordSet = (atom1.cX, atom1.cY, atom2.cX, atom2.cY)
        bond = Bond(coordSet, atom1, atom2)
        data.bondList.append(bond)
    elif data.bondType == 'double':
        coordSet, coordSet2, coordSet3, atom1, atom2 = generateCoordSets(data, 'db')
        bond = DoubleBond(coordSet, coordSet2, coordSet3, atom1, atom2)
        data.bondList.append(bond)
    elif data.bondType == 'triple':
        coordSet, coordSet2, coordSet3, atom1, atom2 = generateCoordSets(data, 'tb')
        bond = TripleBond(coordSet, coordSet2, coordSet3, atom1, atom2)
        data.bondList.append(bond)
        

# =============================================================================
# Represent the connectivity matrix as a dictionary as follows:
# {Atom1 : [(Atom2, 'sb'), (Atom3, 'sb'), (Atom4, 'sb'), (Atom5, 'sb'), 'tetrahedral']...}
# Each atom is a key that connects to other atoms, and we will draw based of the specific geometries
# =============================================================================
def generateConnectivityMatrix(data): #Generates the adjacency matrix for 3D drawing
    connectMatrix = {}
    for atom in data.atoms:
        atomID, tempList = atom.ID, []
        for bond in data.bondList:
            tempAddition = None
            joinedAtomIDs = [atom.ID for atom in bond.connectingAtoms]
            joinedAtomIDSet = set(joinedAtomIDs)
            if atomID in joinedAtomIDSet:
                joinedAtomIDSet.remove(atomID)
                remainingAtom = list(joinedAtomIDSet)[0]
                if isinstance(bond, Bond) and not (isinstance(bond, DoubleBond) or isinstance(bond, TripleBond)):
                    tempAddition = (remainingAtom, 'sb')
                elif isinstance(bond, DoubleBond) and not (isinstance(bond, TripleBond)):
                    tempAddition = (remainingAtom, 'db')
                elif isinstance(bond, TripleBond):
                    tempAddition = (remainingAtom, 'tb')
                tempList.append(tempAddition)
        tempList.append(atom.geometry)
        connectMatrix[atomID] = tempList
    return connectMatrix

#Generate a partial connectivity matrix from the ring data
def genConnectMatrixFromRingData(data):
    partialConnectMatrix = {}
    for i in range(len(data.rings)):
        ring = data.rings[i]
        for j in range(len(ring.connectedAtoms)):
            tempList, currAtom = [], ring.connectedAtoms[j].ID
            prevIndex, nextIndex = j - 1, j + 1
            #Handles indexing out of range
            if nextIndex > len(ring.connectedAtoms) - 1:
                nextIndex = nextIndex % len(ring.connectedAtoms)
            tempList.append((ring.connectedAtoms[prevIndex].ID, 'sb'))
            tempList.append((ring.connectedAtoms[nextIndex].ID, 'sb'))
            tempList.append('ring' + str(i))
            partialConnectMatrix[currAtom] = tempList
    updateValencyInRing(data)
    return partialConnectMatrix

#TODO: Continue here
#Updates the connectivity matrix generated from ring data
def updatePartialConnectivityMatrix(partialRingMatrix, data):
    originalRingMatrixKeys = set(partialRingMatrix.keys())
    for atom in data.atoms:
        atomID = atom.ID
        if atomID in originalRingMatrixKeys:
            currResult = partialRingMatrix[atomID]
            tempLst = []
            for bond in data.bondList:
                joinedAtomIDs = [atm.ID for atm in bond.connectingAtoms]
                joinedAtomIDSet = set(joinedAtomIDs)
                if atomID in joinedAtomIDSet:
                    joinedAtomIDSet.remove(atomID)
                    remainingAtom = list(joinedAtomIDSet)[0]
                    if isinstance(bond, Bond) and not (isinstance(bond, DoubleBond) or isinstance(bond, TripleBond)):
                        tempLst.append((remainingAtom, 'sb'))
                    elif isinstance(bond, DoubleBond) and not (isinstance(bond, TripleBond)):
                        tempLst.append((remainingAtom, 'db'))
                    elif isinstance(bond, TripleBond):
                        tempLst.append((remainingAtom, 'tb'))
            finalResult = currResult + tempLst
            partialRingMatrix[atomID] = finalResult
        elif atomID not in originalRingMatrixKeys:
            tempLst = []
            for bond in data.bondList:
                tempAddition = None
                joinedAtomIDs = [atm.ID for atm in bond.connectingAtoms]
                joinedAtomIDSet = set(joinedAtomIDs)
                if atomID in joinedAtomIDSet:
                    joinedAtomIDSet.remove(atomID)
                    remainingAtom = list(joinedAtomIDSet)[0]
                    if isinstance(bond, Bond) and not (isinstance(bond, DoubleBond) or isinstance(bond, TripleBond)):
                        tempAddition = (remainingAtom, 'sb')
                    elif isinstance(bond, DoubleBond) and not (isinstance(bond, TripleBond)):
                        tempAddition = (remainingAtom, 'db')
                    elif isinstance(bond, TripleBond):
                        tempAddition = (remainingAtom, 'tb')
                    tempLst.append(tempAddition)
            tempLst.append(atom.geometry) #Only add the geometry if the atom is not in a ring
            partialRingMatrix[atomID] = tempLst


def addRingEntries(partialMatrixUpdated, data):
    masterRingList = []
    for ring in data.rings:
        connectedAtoms = ring.connectedAtoms
        idList = [atm.ID for atm in connectedAtoms]
        masterRingList.append(idList)
    partialMatrixUpdated['rings'] = masterRingList


def updateValencyInRing(data): #By default, atoms in a ring will only have, at max, four atoms bonded to it
    for ring in data.rings:
        for atom in ring.connectedAtoms:
            if atom != None:
                atom.numSingleBonds = 2
                updateCurrValency(atom)
    data.ringValencyUpdated = True
    
#Checks if there are solo atoms; molecule fails immediately if this is true
def checkConectivitySoloAtoms(matrix):
    for atom in matrix:
        connectedAtoms = matrix[atom]
        if len(connectedAtoms) == 1:
            return False
    return True

def keyPressed(event, data):
    # use event.char and event.keysym
    if event.char == 'Z':
        if data.cursorMode == 'atomPlacement':
            if len(data.atoms) != 0:
                data.atoms.pop()
        elif data.cursorMode == 'bondPlacement':
            if len(data.bondList) != 0:
                data.bondList.pop()
    elif event.char == 'l':
        data.cursorMode = 'lonePairPlacement'
        data.stringVar.set("Formula: %s    Weight: %0.2f g/mol    Cursor Mode: %s    Atom: %s" % 
                           (data.moleculeFormula, data.moleculeWeight, data.cursorMode, data.atom))
    elif event.char == 'b':
        bondPopup(data)
    elif event.char == 'a':
        data.cursorMode = 'atomPlacement'
        data.stringVar.set("Formula: %s    Weight: %0.2f g/mol    Cursor Mode: %s    Atom: %s" % 
                           (data.moleculeFormula, data.moleculeWeight, data.cursorMode, data.atom))
    elif event.char == 'r':
        partialInit(data)
        updateFormulaAndWeight(data)
    elif event.char == 'i':
        ringPopUp(data)

def redrawAll(canvas, data):
    for bond in data.bondList:
        bond.draw(canvas)
    for ring in data.rings:
        ring.draw(canvas)
    for atom in data.atoms:
        atom.draw(canvas)

#Menu Functions to be called
def resetData(canvas, data):
    partialInit(data)
    canvas.delete(ALL)
    canvas.create_rectangle(0, 0, data.width, data.height,
                            fill='white', width=0)
    redrawAll(canvas, data)
    canvas.update()

#Saving, loading: Pickle, filedialog?
def updateReference(data): #Relinks items together after a pickle file has been loaded
    #Update bond reference first
    for bond in data.bondList:
        connectedAtoms = bond.connectingAtoms
        connectedAtomsIDs = [atom.ID for atom in connectedAtoms]
        for i in range(len(connectedAtomsIDs)):
            currID = connectedAtomsIDs[i]
            for atm in data.atoms:
                if atm.ID == currID:
                    connectedAtoms[i] = atm
        bond.connectingAtoms = connectedAtoms
    #Update ring references next
    for ring in data.rings:
        joinedAtoms = ring.connectedAtoms
        joinedAtomsIDs = [elem.ID for elem in joinedAtoms]
        for j in range(len(joinedAtomsIDs)):
            currentID = joinedAtomsIDs[j]
            for atomos in data.atoms:
                if atomos.ID == currentID:
                    joinedAtoms[j] = atomos
        ring.connectedAtoms = joinedAtoms

def saveLewisStructures(data):
    filePath = os.path.join(os.getcwd(), 'LewisStructs', data.moleculeFormula) + '.p'
    file_out = open(filePath, 'wb')
    pickle.dump(data.atoms, file_out)
    pickle.dump(data.bondList, file_out)
    pickle.dump(data.atomDict, file_out)
    pickle.dump(data.rings, file_out)
    pickle.dump(data.moleculeFormula, file_out)
    pickle.dump(data.moleculeWeight, file_out)
    #The objects are serialized in order, so load them in order
    file_out.close()
    
def saveAsLewisStructures(data):
    tempPop = Tk()
    #Icon taken from https://www.flaticon.com/free-icon/molecule_195742
    tempPop.iconbitmap(os.path.join(os.getcwd(), 'BackgroundImages', 'molecule.ico'))
    tempPop.withdraw()
    initialDir = os.path.join(os.getcwd(), 'LewisStructs')
    tempPop.filename = filedialog.asksaveasfilename(initialdir = initialDir, title = "Save as", filetypes = (("p files","*.p"),("all files","*.*")))
    if tempPop.filename == None or tempPop.filename == '':
        pass
    else:
        file_out = open(tempPop.filename + '.p', 'wb')
        pickle.dump(data.atoms, file_out)
        pickle.dump(data.bondList, file_out)
        pickle.dump(data.atomDict, file_out)
        pickle.dump(data.rings, file_out)
        pickle.dump(data.moleculeFormula, file_out)
        pickle.dump(data.moleculeWeight, file_out) 
        file_out.close()
    tempPop.destroy()

def loadLewisStructures(data, canvas):
    tempPop = Tk()
    #Icon taken from https://www.flaticon.com/free-icon/molecule_195742
    tempPop.iconbitmap(os.path.join(os.getcwd(), 'BackgroundImages', 'molecule.ico'))
    tempPop.withdraw()
    initialDir = os.path.join(os.getcwd(), 'LewisStructs')
    tempPop.filename = filedialog.askopenfile(initialdir = initialDir, title = "Choose your structure", filetypes = (("p files","*.p"),("all files","*.*")))
    if tempPop.filename == None:
        pass
    else:
        partialInit(data)
        path = tempPop.filename.name
        file = open(path, 'rb')
        atomLst = pickle.load(file)
        bondLst = pickle.load(file)
        atomDict = pickle.load(file)
        ringLst = pickle.load(file)
        moleculeFormula = pickle.load(file)
        moleculeWeight = pickle.load(file)
        data.atoms = atomLst
        data.bondList = bondLst
        data.atomDict = atomDict
        data.rings = ringLst
        data.moleculeFormula = moleculeFormula
        data.moleculeWeight = moleculeWeight
        updateReference(data)
        updateFormulaAndWeight(data)
        canvas.delete(ALL)
        canvas.create_rectangle(0, 0, data.width, data.height,
                            fill='white', width=0)
        redrawAll(canvas, data)
        canvas.update()    
    tempPop.destroy()

#Updating cursorMode
def updateCursorMode(flag, data):
    data.cursorMode = flag
    data.stringVar.set("Formula: %s    Weight: %0.2f g/mol    Cursor Mode: %s    Atom: %s" % (data.moleculeFormula, data.moleculeWeight, data.cursorMode, data.atom))


#Updating the atom
def updateAtom(flag, data):
    data.atom = flag
    data.stringVar.set("Formula: %s    Weight: %0.2f g/mol    Cursor Mode: %s    Atom: %s" % (data.moleculeFormula, data.moleculeWeight, data.cursorMode, data.atom))

def checkNoNones(data):
    for ring in data.rings:
        if None in ring.connectedAtoms:
            return False
    return True

#Running ring and non ring visualizations
def ringVisualization(data):
    if len(data.rings) != 0 and checkNoNones(data):
        matrix = genConnectMatrixFromRingData(data)
        updatePartialConnectivityMatrix(matrix, data)
        addRingEntries(matrix, data)
        Geoms.mainFlowWithRings(matrix, data.moleculeFormula, data.moleculeWeight)
    else:
        popupMsg("Error with rings!")

def normVisualization(data):
    if len(data.rings) == 0 and len(data.bondList) != 0:
        matrix = generateConnectivityMatrix(data)
        if checkConectivitySoloAtoms(matrix):
            Geoms.mainFlow(matrix, data.moleculeFormula, data.moleculeWeight)

#Updating Bond Mode
def updateBondMode(flag, data):
    data.bondType = flag

def bondPopup(data):
    popUp = Tk()
    #Icon taken from https://www.flaticon.com/free-icon/molecule_195742
    popUp.iconbitmap(os.path.join(os.getcwd(), 'BackgroundImages', 'molecule.ico'))
    popUp.resizable(width = False, height = False)
    popUp.wm_title('Selecting Bond Type')
    label = Label(popUp, text = 'Please Select the Desired Bond Type', wraplength = 200, font = 'Arial 12')
    label.pack(side = 'top', fill = X, pady = 10)
    B1 = Button(popUp, text = 'Single', command = lambda: [updateBondMode('single', data), updateCursorMode('bondPlacement', data), popUp.destroy()])
    B2 = Button(popUp, text = 'Double', command = lambda: [updateBondMode('double', data), updateCursorMode('bondPlacement', data), popUp.destroy()])
    B3 = Button(popUp, text = 'Triple', command = lambda: [updateBondMode('triple', data), updateCursorMode('bondPlacement', data), popUp.destroy()])
    B1.pack(fill = X), B2.pack(fill = X), B3.pack(fill = X)
    popUp.mainloop()

def isInt(n):
    try:
        int(n)
        return True
    except ValueError:
        return False

def processRingVertexData(text, data):
    if not isInt(text): 
        popupMsg('Illegal Input!')
        return False
    else:
        number = int(text)
        if number > 8: 
            popupMsg('Illegal Input!')
            return False
        elif number < 3:
            popupMsg('Illegal Input!')
            return False
        else:
            data.ringVertices = number
            data.cursorMode = 'ringPlacement'
            data.stringVar.set("Formula: %s    Weight: %0.2f g/mol    Cursor Mode: %s    Atom: %s" % (data.moleculeFormula, data.moleculeWeight, data.cursorMode, data.atom))
            return True

def ringPopUp(data):
    ringPop = Tk()
    #Icon taken from https://www.flaticon.com/free-icon/molecule_195742
    ringPop.iconbitmap(os.path.join(os.getcwd(), 'BackgroundImages', 'molecule.ico'))
    ringPop.resizable(width = False, height = False)
    ringPop.wm_title('Ring Vertices')
    label = Label(ringPop, text = 'Please Enter Number of Atoms (3 < x < 9)', wraplength = 400, font = 'Arial 12')
    label.pack(side = 'top', fill = X, pady = 10)
    entry = Entry(ringPop)
    entry.pack(fill = X)
    B1 = Button(ringPop, text = 'Enter', command = lambda: [ringPop.destroy() if processRingVertexData(entry.get(), data) else None])
    B1.pack(fill = X)
    ringPop.mainloop()

#Gallery popup
def galleryPopUp(data):
    popUp = Tk()
    directoryPath = os.path.join(os.getcwd(), 'MoleculeGallery')
    listOfMolecs = os.listdir(directoryPath)
    popUp.geometry('400x400')
    #Icon taken from https://www.flaticon.com/free-icon/molecule_195742
    popUp.iconbitmap(os.path.join(os.getcwd(), 'BackgroundImages', 'molecule.ico'))
    #Image taken from https://pngio.com/images/png-41028.html
    imagePath = os.path.join(os.getcwd(), 'BackgroundImages', 'Molecule1') + '.png'
    image = Image.open(imagePath)
    photo = ImageTk.PhotoImage(image, master = popUp) #important: set master to the popup
    backgroundLabel = Label(popUp, image = photo)
    backgroundLabel.image = photo #Keeping a reference
    backgroundLabel.place(x = 0, y = 0, relwidth = 1, relheight = 1)
    popUp.resizable(width = False, height = False)
    popUp.wm_title('Gallery Visualizer')
    label = Label(popUp, text = 'Enter a formula with all numbers (e.g. C1H4) and click the visualize button!', wraplength = 400, font = 'Arial 12', bg = None)
    label.pack(fill = X, pady = 10)
    entry = Entry(popUp)
    entry.pack(fill = X)
    B1 = Button(popUp, text = 'Visualize', command = lambda: galleryVisualization(entry.get()))
    B3 = Button(popUp, text = 'Molecule List', command = lambda: moleculeListPopup(data, listOfMolecs))
    B2 = Button(popUp, text = 'Exit', command = lambda: popUp.destroy())
    B1.pack(fill = X), B3.pack(fill = X), B2.pack(fill = X)
    label2 = Label(popUp, text = '''
Some popular choices:
C6N1H11
C2H4
C1H4
C2N1O1H5
C3N1O1H9
C5N1O1H11
C6O1H14
C7H14
C5O1H12
C6N1H15
''', font = 'Arial 10')
    label2.pack(pady = 10)
    popUp.mainloop()


def moleculeListPopup(data, moleculeList):
    popUpPanel = Tk()
    popUpPanel.iconbitmap(os.path.join(os.getcwd(), 'BackgroundImages', 'molecule.ico'))
    popUpPanel.geometry('200x400')
    popUpPanel.resizable(width = False, height = True)
    popUpPanel.wm_title('Molecule List')
    scrollBar = Scrollbar(popUpPanel)
    scrollBar.pack(side = RIGHT, fill = Y)
    listBox = Listbox(popUpPanel, width = 200, height = 400, yscrollcommand = scrollBar.set)
    for i in range(len(moleculeList)):
        listBox.insert(END, moleculeList[i][:-2])
    listBox.pack(side = LEFT, fill = BOTH)
    scrollBar.config(command = listBox.yview)
    popUpPanel.mainloop()



#Taken from: https://www.cs.cmu.edu/~112/notes/notes-recursion.html#permutations
# Problem: given a list, a, find all possible permutations (orderings) of the elements of a
def permutations(a):
    # Base Case: the only permutation of an empty list is the empty list
    if (len(a) == 0):
        return [ [] ]
    else:
        # Recursive Case: remove the first element, then find all possible permutations
        # of the remaining elements. For each permutation, insert a into every possible position.
        partialPermutations = permutations(a[1:])
        allPerms = [ ]
        for subPermutation in partialPermutations:
            for i in range(len(subPermutation) + 1):
                allPerms.append(subPermutation[:i] + [ a[0] ] + subPermutation[i:])
        #returns a 2D list
        return allPerms


#Takes in a 2D list
def combineStringListElements(lst):
    result = ''
    for elem in lst:
        result += elem
    return result

def convertMoleculeFormula(formula):
    result = ''
    for c in formula:
        if c.isalpha():
            result += c.upper()
        else:
            result += c
    return result

def galleryVisualization(moleculeFormula):
    moleculeFormula = convertMoleculeFormula(moleculeFormula)
    directoryPath = os.path.join(os.getcwd(), 'MoleculeGallery')
    setOfMolecules = set(os.listdir(directoryPath))
    moleculeNameSplit = re.split('(\d+)', moleculeFormula)
    moleculeNameSplit.pop() #There will always be an empty string at the end
    atomLst = []
    for i in range(len(moleculeNameSplit)):
        if i % 2 == 0:
            atomAndNum = moleculeNameSplit[i] + moleculeNameSplit[i + 1]
            atomLst.append(atomAndNum)
    allPossibleCombos = permutations(atomLst) #Returns a 2D list
    allPossibleFormulaList = []
    for i in range(len(allPossibleCombos)):
        currList = allPossibleCombos[i]
        result = combineStringListElements(currList)
        allPossibleFormulaList.append(result)
    containedFormula = None
    for formula in allPossibleFormulaList:
        if formula + '.p' in setOfMolecules:
            containedFormula = formula + '.p'
    if containedFormula == None: #Molecule not found
        popupMsg('Illegal Input!')
        return
    else: #molecule found
        completeFilePath = os.path.join(directoryPath, containedFormula)
        bondMatrixFinal, originalMatrix = runGetAnglesAndBondsBase(completeFilePath)
        masterList = []
        for col in range(len(originalMatrix[0])):
            tempList = []
            for row in range(len(originalMatrix)):
                tempList.append(originalMatrix[row][col])
            masterList.append(tempList)
        tempResult = re.split('(\d+)', moleculeFormula)
        tempResult.pop()
        weight = 0
        for i in range(len(tempResult)):
            if i % 2 == 0:
                symbol, num = tempResult[i], tempResult[i + 1]
                elemObj = element(symbol)
                tempWeight = int(num) * elemObj.atomic_weight
                weight += tempWeight
        Geoms.mainFLowWithFinishedData(masterList, bondMatrixFinal, moleculeFormula, weight)

def helpMenuDrawingInfo():
    helpPopUp = Tk()
    helpPopUp.geometry('700x600')
    helpPopUp.wm_title('Instructions for Drawing and General Info')
    #Icon taken from https://www.flaticon.com/free-icon/molecule_195742
    helpPopUp.iconbitmap(os.path.join(os.getcwd(), 'BackgroundImages', 'molecule.ico'))
    helpPopUp.resizable(width = False, height = False)
    #Image taken from http://www.pngmart.com/image/112706
    imagePath = os.path.join(os.getcwd(), 'BackgroundImages', 'Molecule5') + '.png'
    image = Image.open(imagePath)
    image.putalpha(64)
    photo = ImageTk.PhotoImage(image, master = helpPopUp)
    label = Label(helpPopUp, compound = CENTER, text = '''
TO LOAD / SAVE MOLECULES
    1. Go to File > Load to load; you will be prompted to select a file
    2. Go to File > Save to save; saves file to LewisStructs folder, file name is the molecule formula
    3. Go to File > Save as; you will be prompted for a file name
TO PLACE ATOMS
    1. Click anywhere on the canvas to place an atom (default is carbon)
    2. To switch the atom, go to "Atom Options" and select the one you would like
    3. To remove an atom, simply right-click on / near that atom
    4. To go back to atom placement mode, do Cursor > Atom Placement, or press 'a'
    5. If you try to add more bonds to an atom than allowed, you will be told "Maximum Valency Exceeded!"
TO CONNECT ATOMS
    1. Change to bond mode by doing Cursor > Bond Placement, or by pressing 'b'
    2. You will be prompted to enter the bond type (single, double, triple)
    3. Click on the first atom. It will turn green
    4. Click on the second atom, and the bond will form
    5. To de-select atoms, left click anywhere on the empty canvas
    6. To change bond order, press 'b' or go to Cursor > Bond Placement and select the appropriate option
    7. To remove bonds, right-click near the bond
TO MAKE RINGS
    DO NOT make rings with normal bonds! Instantiate rings through ring placement mode only.
    1. Go to Cursor > Ring Placement, or press 'i'
    2. You will be prompted to enter the number of atoms (3 < x < 8)
    3. If you enter anything outside that range, or a non-integer, you will be told "Illegal Input!"
    4. Once you place the ring, you can populate the vertices with atoms by clicking in atom placement mode
    5. Make sure that all vertices are populated before trying to visualize!
    6. To remove rings, right-click near the center of the ring
TO ADD LONE PAIRS
    1. Go to Cursor > Lone Pair Placement, or press 'l'
    2. Click on the atom you would like to place lone pairs on
    3. If you go overboard, you will be told "Maximum Valency Exceeded!"
TO RESET CANVAS
    1. Go to File > Reset, or press 'r'
TO EXIT THE APPLICATION
    1. Go to File > Exit, or simply close the main window
    2. Application windows should be closed individually
    3. The browser tab should be closed last
''', font = 'Arial 10 bold', fg = 'black',  image = photo,  justify = LEFT)
    label.pack()
    helpPopUp.mainloop()


def helpMenu3DInfo():
    helppop = Tk()
    helppop.geometry('550x350')
    helppop.wm_title('Instructions for 3D Navigation')
    #Icon taken from https://www.flaticon.com/free-icon/molecule_195742
    helppop.iconbitmap(os.path.join(os.getcwd(), 'BackgroundImages', 'molecule.ico'))
    helppop.resizable(width = False, height = False)
    #Image taken from http://www.scienceclarified.com/Ma-Mu/Molecule.html
    imagePath = os.path.join(os.getcwd(), 'BackgroundImages', 'molecule-3000') + '.jpg'
    image = Image.open(imagePath)
    image.putalpha(64)
    photo = ImageTk.PhotoImage(image, master = helppop)
    label = Label(helppop, compound = CENTER, text = '''
VISUALIZING THE MOLECULE
    1. Make sure that when you visualize the molecule, you select the right option
    2. If the molecule contains rings, select the option that allows for rings
    3. If no rings, select the appropriate option that says no rings
TO ROTATE / ZOOM THE MOLECULE
    1. Right-click on the scene, hold down, and drag
    2. To zoom in and out, use the mouse wheel or left-click and drag
TO CLEAR THE SCENE
    1. Go to Visualize! > Clear Screen
    2. The screen also clears between visualizations of molecules
CAPTURING SCREENSHOTS
    1. Click on the "Screenshot!" button in the browser window
    2. All screenshots are sent to your downloads folder
CHECKING COLOR
    1. Go to Help > Element Color Code to pull up a color code chart
    2. Drag and expand the color code chart to move it more easily
EXITING THE LOCALHOST
    1. Close the application windows individually first
    2. Close the browser tab last
''', font = 'Arial 10 bold', fg = 'black', image = photo, justify = LEFT)
    label.pack()
    helppop.mainloop()

def helpMenuGalleryVisuals():
    helpMenuPop = Tk()
    helpMenuPop.geometry('550x200')
    #Icon taken from https://www.flaticon.com/free-icon/molecule_195742
    helpMenuPop.iconbitmap(os.path.join(os.getcwd(), 'BackgroundImages', 'molecule.ico'))
    helpMenuPop.wm_title('Instructions for Gallery Visualizations')
    helpMenuPop.resizable(width = False, height = False)
    #Image taken from http://www.scienceclarified.com/Ma-Mu/Molecule.html
    imagePath = os.path.join(os.getcwd(), 'BackgroundImages', 'molecule-3000') + '.jpg'
    image = Image.open(imagePath)
    image.putalpha(64)
    photo = ImageTk.PhotoImage(image, master = helpMenuPop)
    label = Label(helpMenuPop, compound = CENTER, text = '''
VISUALIZING USING THE GALLERY
    IMPORTANT: Please wait for one rendering to finish before requesting another
    one! Spamming might lead to program freezes
    1. Click on the Gallery Visualizer button on the menu
    2. Enter a chemical formula
    3. The elements can be randomly ordered, and it will still work (C1H4 == H4C1)
    4. Click Visualize to see the molecule rendered in 3D space
    5. Click Exit to close the visualization panel
    6. You can view all possible formulas by clicking on the Molecule List button
''', font = 'Arial 10 bold', fg = 'black', image = photo, justify = LEFT)
    label.pack()
    helpMenuPop.mainloop()
    pass

def colorCodeChart():
    colorCodeChart = Tk()
    colorCodeChart.geometry('300x400')
    colorCodeChart.resizable(width = False, height = False)
    colorCodeChart.wm_title('Color Code for Atoms')
    #Icon taken from https://www.flaticon.com/free-icon/molecule_195742
    colorCodeChart.iconbitmap(os.path.join(os.getcwd(), 'BackgroundImages', 'molecule.ico'))
    colorCodeChart.resizable(width = False, height = False)
    imagePath = os.path.join(os.getcwd(), 'BackgroundImages', 'colorCode.png')
    image = Image.open(imagePath)
    photo = ImageTk.PhotoImage(image, master = colorCodeChart)
    backgroundLabel = Label(colorCodeChart, image = photo)
    backgroundLabel.image = photo
    backgroundLabel.pack()
    colorCodeChart.mainloop()

def FAQs():
    faqPanel = Tk()
    faqPanel.wm_title('Frequently Asked Questions')
    faqPanel.geometry('600x300')
    #Icon taken from https://www.flaticon.com/free-icon/molecule_195742
    faqPanel.iconbitmap(os.path.join(os.getcwd(), 'BackgroundImages', 'molecule.ico'))
    faqPanel.resizable(width = False, height = False)
    #Image taken from https://pngio.com/images/png-41028.html
    imagePath = os.path.join(os.getcwd(), 'BackgroundImages', 'Molecule1.png')
    image = Image.open(imagePath)
    image.putalpha(64)
    photo = ImageTk.PhotoImage(image, master = faqPanel)
    label = Label(faqPanel, compound = CENTER, text = '''
Q: Why do some of the molecules in the gallery have disjoints?
    A: Bond data are imperfect, so some bond lengths may have been missed.
Q: Why did you decide to do this for your project?
    A: I'm a chemistry major. I know nothing else.
Q: How many hours did you REALLY put into this project?
    A: If I'm being honest, I lost track after TP2. I think around 40?
Q: Are you proud of your work?
    A: Considering what I put into it, I think the end product is decent
Q: What was the hardest part?
    A: Figuring out all the nuances of the 2D Lewis structures
Q: What was the easiest part?
    A: The UI, relative to everything else, was a lot easier
Q: Did you learn a lot while doing this project?
    A: Yes.
Q: What did it cost?
    A: Everything.
''', font = 'Arial 10 bold', fg = 'black', image = photo, justify = LEFT)
    label.pack()
    faqPanel.mainloop()

def importantChemicalInformation():
    chemInfo = Tk()
    chemInfo.wm_title('Important Chemical Information')
    chemInfo.geometry('600x200')
    #Icon taken from https://www.flaticon.com/free-icon/molecule_195742
    chemInfo.iconbitmap(os.path.join(os.getcwd(), 'BackgroundImages', 'molecule.ico'))
    chemInfo.resizable(width = False, height = False)
    #Image taken from http://www.scienceclarified.com/Ma-Mu/Molecule.html
    imagePath = os.path.join(os.getcwd(), 'BackgroundImages', 'molecule-3000.jpg')
    image = Image.open(imagePath)
    image.putalpha(64)
    photo = ImageTk.PhotoImage(image, master = chemInfo)
    label = Label(chemInfo, compound = CENTER, text = '''
BONDING INFORMATION
    
    Boron (B):     3 bonds max
    Carbon (C):    4 bonds max
    Nitrogen (N):  3 bonds max
    Oxygen (O):    2 bonds max
    Fluorine (F):  1 bond max
    
    Elements in period 2 can have up to 4 lone pairs. Elements starting in period 
    3 and onward can form up to 6 bonds or have 6 lone pairs, with the maximum
    geometry of octahedral. Double and triple bonds count as 2 and 3 bonds, respectively.
''', font = 'Arial 10 bold', fg = 'black', image = photo, justify = LEFT)
    label.pack()
    chemInfo.mainloop()
    
    
####################################
# use the run function as-is
####################################
#Code taken + modified from https://www.cs.cmu.edu/~112/notes/notes-animations-part1.html
def run(width=300, height=300):
    def redrawAllWrapper(canvas, data):
        canvas.delete(ALL)
        canvas.create_rectangle(0, 0, data.width, data.height,
                                fill='white', width=0)
        redrawAll(canvas, data)
        canvas.update()

    def mousePressedWrapper(event, canvas, data):
        mousePressed(event, data)
        redrawAllWrapper(canvas, data)

    def keyPressedWrapper(event, canvas, data):
        keyPressed(event, data)
        redrawAllWrapper(canvas, data)
    
    def mouseRightClickWrapper(event, canvas, data):
        mouseRightClick(event, data)
        redrawAllWrapper(canvas, data)

    # Set up data and call init
    class Struct(object): pass
    data = Struct()
    data.width = width
    data.height = height
    root = Tk()
    root.resizable(width=False, height=False) # prevents resizing window
    root.title('MoleCool')
    #Icon taken from https://www.flaticon.com/free-icon/molecule_195742
    root.iconbitmap(os.path.join(os.getcwd(), 'BackgroundImages', 'molecule.ico'))
    init(data)
    #Bottom Frame Tool bar
    botFrame = Frame(root)
    botFrame.pack(side = BOTTOM)
    #Formula and weight labels
    label = Label(botFrame, textvariable = data.stringVar, anchor = W)
    label.pack()
    # create the root and the canvas
    canvas = Canvas(root, width=data.width, height=data.height)
    canvas.configure(bd=0, highlightthickness=0)
    canvas.pack()
    #Setup Menu
    #Wrap everything in lambdas!
    menu = Menu(root)
    root.config(menu = menu)
    
    fileMenu = Menu(menu, tearoff = 0)
    menu.add_cascade(label = 'File', menu = fileMenu)
    fileMenu.add_command(label = 'Save', command = lambda: saveLewisStructures(data))
    fileMenu.add_command(label = 'Save As', command = lambda: saveAsLewisStructures(data))
    fileMenu.add_command(label = 'Load', command = lambda: loadLewisStructures(data, canvas))
    fileMenu.add_separator()
    fileMenu.add_command(label = 'Reset (r)', command = lambda: resetData(canvas, data)) 
    fileMenu.add_command(label = 'Exit', command = lambda: root.destroy()) 
    
    cursorMenu = Menu(menu, tearoff = 0)
    menu.add_cascade(label = 'Cursor', menu = cursorMenu)
    cursorMenu.add_command(label = 'Bond Placement (b)', command = lambda: bondPopup(data))
    cursorMenu.add_command(label = 'Atom Placement (a)', command = lambda: updateCursorMode('atomPlacement', data))
    cursorMenu.add_command(label = 'Ring Placement (i)', command = lambda: ringPopUp(data))
    cursorMenu.add_command(label = 'Lone Pair Placement (l)', command = lambda: updateCursorMode('lonePairPlacement', data))

    
    atomMenu = Menu(menu, tearoff = 0)
    menu.add_cascade(label = 'Atom Options', menu = atomMenu)
    atomMenu.add_command(label = 'Hydrogen (H)', command = lambda: updateAtom('H', data)) 
    atomMenu.add_command(label = 'Fluorine (F)', command = lambda: updateAtom('F', data)) 
    atomMenu.add_command(label = 'Chlorine (Cl)', command = lambda: updateAtom('Cl', data))
    atomMenu.add_command(label = 'Bromine (Br)', command = lambda: updateAtom('Br', data))
    atomMenu.add_command(label = 'Iodine (I)', command = lambda: updateAtom('I', data))
    atomMenu.add_command(label = 'Oxygen (O)', command = lambda: updateAtom('O', data))
    atomMenu.add_command(label = 'Sulfur (S)', command = lambda: updateAtom('S', data))
    atomMenu.add_command(label = 'Selenium (Se)', command = lambda: updateAtom('Se', data))
    atomMenu.add_command(label = 'Tellurium (Te)', command = lambda: updateAtom('Te', data))
    atomMenu.add_command(label = 'Nitrogen (N)', command = lambda: updateAtom('N', data))
    atomMenu.add_command(label = 'Phosphorus (P)', command = lambda: updateAtom('P', data))
    atomMenu.add_command(label = 'Aresenic (As)', command = lambda: updateAtom('As', data))
    atomMenu.add_command(label = 'Antimony (Sb)', command = lambda: updateAtom('Sb', data))
    atomMenu.add_command(label = 'Carbon (C)', command = lambda: updateAtom('C', data))
    atomMenu.add_command(label = 'Silicon (Si)', command = lambda: updateAtom('Si', data))
    atomMenu.add_command(label = 'Germanium (Ge)', command = lambda: updateAtom('Ge', data))
    atomMenu.add_command(label = 'Boron (B)', command = lambda: updateAtom('B', data))
    atomMenu.add_command(label = 'Xenon (Xe)', command = lambda: updateAtom('Xe', data))
    atomMenu.add_command(label = 'Krypton (Kr)', command = lambda: updateAtom('Kr', data))
    
    visualizationMenu = Menu(menu, tearoff = 0)
    menu.add_cascade(label = 'Visualize!', menu = visualizationMenu)
    visualizationMenu.add_command(label = 'Visualize Molecule, Contains Rings', command = lambda: ringVisualization(data))
    visualizationMenu.add_command(label = 'Visualize Molecule, Does Not Contain Rings', command = lambda: normVisualization(data))
    visualizationMenu.add_command(label = 'Clear Screen', command = lambda: Geoms.clearScreen('', 0))

    menu.add_command(label = 'Gallery Visualizer', command = lambda: galleryPopUp(data))
    
    helpMenu = Menu(menu, tearoff = 0)
    menu.add_cascade(label = 'Help', menu = helpMenu)
    helpMenu.add_command(label = 'Element Color Code', command = lambda: colorCodeChart())
    helpMenu.add_command(label = 'Important Chemical Information', command = lambda: importantChemicalInformation())
    helpMenu.add_command(label = 'Drawing and General Info', command = lambda: helpMenuDrawingInfo())
    helpMenu.add_command(label = '3D Visuals Info', command = lambda: helpMenu3DInfo())
    helpMenu.add_command(label = 'Gallery Browsing Info', command = lambda: helpMenuGalleryVisuals())
    helpMenu.add_command(label = 'FAQs', command = lambda: FAQs())
    
    # set up events
    canvas.bind("<Button-1>", lambda event:
                            mousePressedWrapper(event, canvas, data))
    root.bind("<Key>", lambda event:
                            keyPressedWrapper(event, canvas, data))
    canvas.bind("<Button-3>", lambda event: mouseRightClickWrapper(event, canvas, data))
    redrawAllWrapper(canvas, data)
    # and launch the app
    root.mainloop()  # blocks until window is closed
    print("bye!")

