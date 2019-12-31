# -*- coding: utf-8 -*-
#File for calculating 3D geometries from data taken from qm7
import math
import numpy as np
import xlrd
from copy import deepcopy
import re
import os
import pickle

# =============================================================================
# Helper Functions / variables

def bondLengthDistance(x1,y1,z1,x2,y2,z2):
    return (((x2 - x1)**2) + ((y2-y1)**2) + ((z2-z1)**2))**0.5

BohrRadius = 0.529177 #Given in angstroms

def angleCalculationbyCosineLaw(a,b,c):
    return (((a**2) + (b**2)) - (c**2)) / (2*a*b)

def matrixConversion(excelmat):
    result = np.empty((excelmat.nrows, excelmat.ncols),dtype = object)
    for row in range(excelmat.nrows):
        for col in range(excelmat.ncols):
            result[row,col] = excelmat.cell_value(row,col)
    return result

def createExcelMat(path):
    temporaryMat = xlrd.open_workbook(path)
    result = temporaryMat.sheet_by_index(0)
    return result

def dimensionCorrection(originalMatrix):
     multiplier = len(originalMatrix[0])
     if len(originalMatrix) == 4: return originalMatrix
     elif len(originalMatrix) == 2:
         originalMatrix = np.append(originalMatrix, [[0] * multiplier], axis = 0)
         originalMatrix = np.append(originalMatrix, [[0] * multiplier], axis = 0)
     elif len(originalMatrix) == 3:
         originalMatrix = np.append(originalMatrix, [[0] * multiplier], axis = 0)
     return originalMatrix
             
# =============================================================================

class moleculeAndVariables: #Object that initalizes the constants for this molecule
    def __init__(self, path):
        self.originalMatrixPath = path
        objFile = open(self.originalMatrixPath, 'rb')
        matrix = pickle.load(objFile)
        self.unconvertedOriginalMatrix = list(matrix)
        self.originalMatrix = moleculeAndVariables.checkAndFix(self.unconvertedOriginalMatrix)
        self.unconvertedOriginalMatrix = deepcopy(self.originalMatrix)
        self.originalMatrix = np.array(self.originalMatrix, dtype = object)
        #Data for bond lengths taken from http://www.wiredchemist.com/chemistry/data/bond_energies_lengths.html
        self.cutoff = os.getcwd() + os.sep + 'TabulationsofCutoffBondLengths2.xlsx'
        self.originalMatrix = dimensionCorrection(self.originalMatrix)
        bondLengthCutoffMatrix = createExcelMat(self.cutoff)
        self.workingBondLengthCutoffMatrix = matrixConversion(bondLengthCutoffMatrix)
        self.atomCount = len(self.originalMatrix[0])
        self.numRows = 4
    
    @staticmethod
    def checkAndFix(matrix):
        result = deepcopy(matrix)
        for i in range(len(result)):
            if result[i] == []:
                result[i] = [0] * len(result[0])
        return result
                
def initializeMatrices(path): #Called by main func to initalize the matrices
    matrices = moleculeAndVariables(path)
    return matrices

def createStringToNumDict(matrix): #Takes the original matrix and creates the string to num dictionary
    stringToNumDict = dict()
    for i in range(len(matrix[0])):
        stringToNumDict[matrix[0,i]] = float(i)
    return stringToNumDict

def createBondLengthCutoffDict(cutoffMat): #Creates a dictionary of the bondlengthcutoffs
    bondLengthCutoffs = dict() 
    for row in range(1, len(cutoffMat)):
        key , value = cutoffMat[row,0] , cutoffMat[row,1:]
        bondLengthCutoffs[key] = value
    return bondLengthCutoffs


def createEmptyBondMatrix(matrix): #Creates the empty bond matrix of specified number of rows
    workingMatrix = matrix
    workingMatrixDimensions = workingMatrix.shape
    headingList = workingMatrix[0,:].tolist()
    atomString = ''
    for i in range(len(headingList)):
        atomString += headingList[i]
    numberHydrogens = atomString.count('H') #Gives the # of Hydorgens in the molecule, which necessarily cannot form bonds with each other
    numberNonHydrogens = workingMatrixDimensions[1] - numberHydrogens #workingMatrixDimensions[1] is equivalent to the total number of atoms
    bondLengthMatrixRows = 0
    #You start with the maximum number of bonds (the first atom can bond to all atoms after it) and work down by one each time. For example,
    #   in ethane, you have C2H6. The first carbon can form 7, and the second carbon can form 6 (we already counted C1C2!), and we continue
    #   adding by one less until we finish with the number of non-hydrogen atoms
    for i in range(workingMatrixDimensions[1] - 1, (workingMatrixDimensions[1]-numberNonHydrogens)-1, -1): 
        bondLengthMatrixRows += i
    bondLengthMatrix = np.empty((bondLengthMatrixRows, 2), dtype = object)
    return bondLengthMatrix, workingMatrixDimensions


def addInBondLengths(bondLengthMatrix, dimensions, funcMatrix): #Destructively adds things into the bond length function
    place = 0
    for i in range(dimensions[1]):
        if 'H' in funcMatrix[0,i]:
            pass
        else:
            for j in range(i + 1, dimensions[1]):
                index = funcMatrix[0,i] + funcMatrix[0,j]
                bondLengthMatrix[place,0] = index
                X1,X2 = funcMatrix[1,i], funcMatrix[1,j]
                Y1,Y2 = funcMatrix[2,i], funcMatrix[2,j]
                Z1,Z2 = funcMatrix[3,i], funcMatrix[3,j] 
                BondLength = bondLengthDistance(X1,Y1,Z1,X2,Y2,Z2)
                bondLengthMatrix[place,1] = BondLength * BohrRadius
                place += 1
    pass

def removeBadRows(bondLengthMatrix, bondLengthCutoffs): #THIS IS THE BLOCK THAT SORTS THROUGH BOND LENGTHS
    badRows, numberBadRows = [], 0
    for row in range(len(bondLengthMatrix)): #For all bonds that exceed the bond length, the value is set to None
        Checkmatrix = []
        for character in bondLengthMatrix[row,0]:
            if character.isalpha():
                Checkmatrix.append(character)
        Checkstring = Checkmatrix[0] + Checkmatrix[1]
        Threshold = bondLengthCutoffs[Checkstring]
        if (math.isclose(bondLengthMatrix[row, 1], Threshold[0], abs_tol = 1e-1) == False) and\
        (math.isclose(bondLengthMatrix[row, 1], Threshold[1], abs_tol = 1e-1) == False) and\
        (math.isclose(bondLengthMatrix[row, 1], Threshold[2], abs_tol = 1e-1) == False):
            badRows.append(row)
            numberBadRows += 1
    numberGoodRows = bondLengthMatrix.shape[0] - numberBadRows #Takes the number of total rows and subtracts out the number of bad rows from the matrix
    bondLengthMatrixFinal = np.empty((numberGoodRows, 2), dtype = object)
    rowCount = 0

    for row in range(len(bondLengthMatrix)): #Brings all legal bonds into a final matrix
        if row not in badRows:
            bondLengthMatrixFinal[rowCount, 0], bondLengthMatrixFinal[rowCount, 1] = bondLengthMatrix[row, 0], bondLengthMatrix[row, 1]
            checkMatrix = []
            for character in bondLengthMatrix[row, 0]:
                if character.isalpha():
                    checkMatrix.append(character)
            checkString = checkMatrix[0] + checkMatrix[1]
            Threshold = bondLengthCutoffs[checkString]
            if math.isclose(bondLengthMatrix[row, 1], Threshold[0], abs_tol = 1e-1):
                bondLengthMatrixFinal[rowCount][0] = bondLengthMatrixFinal[rowCount][0] + 'sb'
            elif math.isclose(bondLengthMatrix[row, 1], Threshold[1], abs_tol = 1e-1):
                bondLengthMatrixFinal[rowCount][0] = bondLengthMatrixFinal[rowCount][0] + 'db'
            elif math.isclose(bondLengthMatrix[row, 1], Threshold[2], abs_tol = 1e-1):
                bondLengthMatrixFinal[rowCount][0] = bondLengthMatrixFinal[rowCount][0] + 'tb'
            rowCount += 1
    return bondLengthMatrixFinal
        

def getAngleAndBondsBase(matrices, stringToNumDict):
    funcMatrix = np.array(matrices.originalMatrix, dtype = object)
    for i in range(len(funcMatrix[0])): 
        for key in stringToNumDict:
            if stringToNumDict[key] == float(i):
                funcMatrix[0,i] = key
    bondLengthCutoffs = createBondLengthCutoffDict(matrices.workingBondLengthCutoffMatrix)
    bondMatrixReturns = createEmptyBondMatrix(funcMatrix)
    bondMatrix, workingMatrixDimensions = bondMatrixReturns[0], bondMatrixReturns[1]
    addInBondLengths(bondMatrix, workingMatrixDimensions, funcMatrix)
    bondMatrixFinal = removeBadRows(bondMatrix, bondLengthCutoffs)
    bondLengthTemplate = bondMatrixFinal[:, 0].tolist()
    return bondMatrixFinal, bondLengthTemplate


def runGetAnglesAndBondsBase(path):
    matrices = initializeMatrices(path)
    stringToNumDict = createStringToNumDict(matrices.originalMatrix)
    matrices.stringToNumDict = stringToNumDict
    output = getAngleAndBondsBase(matrices, stringToNumDict)
    bondMatrixFinal, bondLengthTemplate = output[0], output[1]
    return bondMatrixFinal, matrices.unconvertedOriginalMatrix
    
