'''
Description: This script will randomly alter and "mistune" a supplied .bdf
file, altering both mass and stiffness properties. It uses pyNastran to read 
.bdf file info, alter the properties, and write a the altered file as a new .bdf

Author: Braden Frigoletto

Date: Sep. 2020
'''

# Python imports
import sys
import os
import numpy as np
import pyNastran 
from pyNastran.bdf.bdf import BDF

#------------------------------------------------------------------------------
# function for mistuning BDF models
#------------------------------------------------------------------------------
def mistuneBDF(bdfName, numChanges, changeType, altInt):
	'''
	Inputs:
	bdfName - (string) file name of the .bdf file that will be mistuned
	numChanges - (int) total number of fields that will be altered in the .bdf
	             Must be <= the total number of fields in the model
	changeType - (string) "stiffness", "mass", or "both", selects the kinds of 
	             model properties that will be altered
	altInt - (int) value of the mistune ID for the file e.g. "test_Alt_7.bdf"

	Example: 
		originalBDF = "test.bdf"
		changesToMake = 6
		changeType = "both"
		for i in range(3):
			mistuneBDF(originalBDF, changesToMake, changeType, i+1)

		This creates an output directory named "mistune_output" which is 
		populated with "test_Alt_1.bdf", "test_Alt_2.bdf", and "test_Alt_3.bdf"

		Each alt file would contain a total of 6 changes to either a stiffness
		property or a mass property. Each alt file is mistuned differently, and 
		the number of stiffness and mass alterations may vary between them. 

		To control the level of the alteration, adjust the variable 'factor'.
		Note that this variable affects different properties in different ways.
		E.g. factor = 0.1 can yield up to a 10% change in a material density 
		value, but can yield up to a 46% change in a beam bending stiffness. 
		This is because the factor is applied to the beam cross-section 
		dimensions and not the stiffness value itself.

		You can also control which properties are more likely to be changed
		by changing the 'threshold' value in each of the following functions to 
		suit your needs. A threshold value of 50 means that both outcomes have
		equal an chance of being selected.
	'''

	# error check on inputs
	if bdfName[-4::] != '.bdf':
		raise NameError("%s is not a .bdf file" % bdfName)

	if numChanges <= 0:
		raise ValueError("Must input a positive integer value for number of changes to make")

	if changeType not in ["stiffness", "mass", "both"]:
		raise ValueError("Type of change to make %s is not \"stiffness\", \"mass\", or \"both\"" % changeType)
	
	
	# create an output directory for the results if one doesn't exist
	outputDir = "mistune_output"
	if not os.path.isdir(outputDir):
		os.mkdir(outputDir)

	# open the original .bdf file
	model = BDF()
	model.read_bdf(bdfName, punch=True)
	print(model.get_bdf_stats())

	# get the number of properties of each type - separate mass and stiffness
	elemPropKeys = list(model.properties.keys())
	matStiffPropKeys = list(model.materials.keys())
	matMassPropKeys = list(model.materials.keys())
	conMassKeys = list(model.masses.keys())

	numElemProps = len(elemPropKeys)
	numMatStiffProps = len(matStiffPropKeys)
	numMatMassProps = len(matMassPropKeys)
	numConMasses = len(conMassKeys)

	# check the total number of properties available to change
	totChangesPossible = 0
	if changeType == "stiffness":
		totChangesPossible += (numElemProps + numMatStiffProps)
	elif changeType == "mass":
		totChangesPossible += (numConMasses + numMatMassProps)
	elif changeType == "both":
		totChangesPossible += (numElemProps + numConMasses + numMatStiffProps + numMatMassProps)

	# raise error if user tries to change more than what is possible
	if numChanges > totChangesPossible:
		raise ValueError("Desired number of changes (%i) is larger than %i available to change" % (numChanges, totChangesPossible))
	
	# how to alter info in a CONM2
	# specified by the EID number
	# print(model.masses[366].X) # access to the offset from grid point
	# print(model.masses[366].mass) # access to the mass value
	# print(model.masses[366].I) # access to the inertia values

	# TO DO: customize the factor for each type of change
	# set boundary for the magnitude of the changes that will be made
	factor = 0.05 

	# run through loop and alter the .bdf entries
	for i in range(numChanges):

		if changeType == "stiffness":
			alterStiffness(model, factor, elemPropKeys, matStiffPropKeys)

		elif changeType == "mass":
			alterMass(model, factor, conMassKeys, matMassPropKeys)

		elif changeType == "both":
			alterBoth(model, factor, elemPropKeys, conMassKeys, matStiffPropKeys, matMassPropKeys)

	# write the changed model to a new file in the output directory
	if bdfName.endswith('.bdf'):
		bdfName = bdfName[:-4]
	bdfFilenameOut = outputDir + '/' + bdfName + '_Alt_' + str(altInt) + '.bdf'
	model.write_bdf(bdfFilenameOut)

	return

#------------------------------------------------------------------------------
# function for altering a stiffness value
#------------------------------------------------------------------------------
def alterStiffness(model, factor, elemPropKeys, matStiffPropKeys):

	# set probability threshold between random options 
	threshold = 50

	# determine a random bounded percentage to alter a property
	percent = np.random.uniform(-1,1)
	changeFactor = 1 + (percent * factor)
	
	# count available changes for each type
	numElemPropLeft = len(elemPropKeys) 
	numMatStiffPropLeft = len(matStiffPropKeys)

	# use random number to decide whether to alter the element property or the 
	# material property to change stiffness
	cont = True
	while cont:
		val = np.random.uniform(0,100)
		if val < threshold and numElemPropLeft > 0:
			# alter the element property (DIMENSIONS OF PBEAML ONLY FOR NOW!!!)

			# select a random elem property in that set
			key_index = -1
			while key_index not in elemPropKeys:
				key_index = np.random.randint(min(elemPropKeys), max(elemPropKeys)+1)
			print("altering element property [%i] dimensions by a factor of: %f" % (key_index, changeFactor))
			
			# get the dimensions of the beam property and multiply all of them by the factor
			dimensions = model.properties[key_index].dim
			model.properties[key_index].dim = [changeFactor * dimensions[i] for i in range(len(dimensions))]
			
			# remove that elem property from the list of available properties
			elemPropKeys.remove(key_index)
			
			cont = False
		
		elif val >= threshold and numMatStiffPropLeft > 0:
			# alter the material property (E,G,nu of MAT1 ONLY FOR NOW!!!)

			# select a random mat property in that set
			key_index = -1
			while key_index not in matStiffPropKeys:
				key_index = np.random.randint(min(matStiffPropKeys), max(matStiffPropKeys)+1)
			print("altering material [%i] stiffness properties by a factor of: %f" % (key_index, changeFactor))
			
			# get the material stiffness properties and multiply all of them by the factor
			# Note that the G calulcation seems to be unecessary as it is not used in the resulting file
			altE = changeFactor * model.materials[key_index].e
			altNu = changeFactor * model.materials[key_index].nu
			_, altG, _ = model.materials[key_index].set_E_G_nu(altE, None, altNu)
			model.materials[key_index].e = altE
			model.materials[key_index].g = altG
			model.materials[key_index].nu = altNu
			
			# remove that stiffness material property from the list of available properties
			matStiffPropKeys.remove(key_index)

			cont = False

	return 

#------------------------------------------------------------------------------
# function for altering a mass property
#------------------------------------------------------------------------------
def alterMass(model, factor, conMassKeys, matMassPropKeys):

	# set probability threshold between random options 
	threshold = 50

	# determine a random bounded percentage to alter a property
	percent = np.random.uniform(-1,1)
	changeFactor = 1 + (percent * factor)

	# count available changes for each type
	numconMassLeft = len(conMassKeys) 
	numMatMassPropLeft = len(matMassPropKeys)

	# use random number to decide whether to alter the con mass property or the 
	# material property to change mass
	cont = True
	while cont:
		val = np.random.uniform(0,100)
		if val < threshold and numconMassLeft > 0:
			# alter the con mass property (mass, inertia, or offset - CONM2)

			# select a random con mass in that set
			key_index = -1
			while key_index not in conMassKeys:
				key_index = np.random.randint(min(conMassKeys), max(conMassKeys)+1)
			print("altering the mass value of con mass [%i] by a factor of: %f" % (key_index, changeFactor))
			
			# alter the mass of the con mass (can later add functionality to alter offset distances or inertias)
			model.masses[key_index].mass *= changeFactor

			# remove that con mass from the list of available con masses
			conMassKeys.remove(key_index)

			cont = False

		elif val >= threshold and numMatMassPropLeft > 0:
			# alter the material property (density - MAT1)

			# select a random material mass property in that set
			key_index = -1
			while key_index not in matMassPropKeys:
				key_index = np.random.randint(min(matMassPropKeys), max(matMassPropKeys)+1)
			print("altering material [%i] density property by a factor of: %f" % (key_index, changeFactor))
			
			# alter the density value for the material property
			model.materials[key_index].rho *= changeFactor
			
			# remove that material mass property from the list of available properties
			matMassPropKeys.remove(key_index)
			
			cont = False
		
	return

#------------------------------------------------------------------------------
# function for altering both mass and stiffness properties
#------------------------------------------------------------------------------
def alterBoth(model, factor, elemPropKeys, conMassKeys, matStiffPropKeys, matMassPropKeys):

	# set probability threshold between random options 
	threshold = 50

	# count available changes for each type
	numStiffnessLeft = len(elemPropKeys) + len(matStiffPropKeys)
	numMassLeft = len(conMassKeys) + len(matMassPropKeys)

	# use random number to decide whether to alter the stiffness or mass 
	cont = True
	while cont:
		val = np.random.uniform(0,100)

		if val < threshold and numStiffnessLeft > 0:
			# alter the stiffness
			alterStiffness(model, factor, elemPropKeys, matStiffPropKeys)
			cont = False

		elif val >= threshold and numMassLeft > 0:
			# alter the mass
			alterMass(model, factor, conMassKeys, matMassPropKeys)
			cont = False

	return
#------------------------------------------------------------------------------

# Test
test_bdf = "beam_model.bdf"
for i in range(5):
	mistuneBDF(test_bdf,15,"mass",i+1)
