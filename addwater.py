import numpy as np
import MDAnalysis
import MDAnalysis.analysis
import sys
import scipy.cluster

######################################################################################################################################################################
######################################################################################################################################################################
def unitvector(v):
	normal = np.linalg.norm(v)
	UV = v/normal
	return UV
######################################################################################################################################################################
######################################################################################################################################################################
def writewaterfile(filename, watercoods):

	numwater = watercoods.shape[0]
	f1 = open(filename,'w')

	for j in xrange(0,numwater):
		header = 'HETATM'
		serial = j+1
		name = 'OW'
		resname = 'SOL'
		chainID = 'A'
		resSeq = j+1
		icode = ' '
		occupancy = 1.0
		tempfactor = 0.0
		x = watercoods[j,0]
		y = watercoods[j,1]
		z = watercoods[j,2]

		f1.write("%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n" %(header,serial,name,icode,resname,chainID,resSeq,icode,x,y,z,occupancy,tempfactor))

	f1.close()
######################################################################################################################################################################
######################################################################################################################################################################
def carbonylorcarboxyl(allligand,index,bond_dist):

	allligandcoods = allligand.positions
	ocoods = np.zeros((1,3), dtype = float)
	ocoods[0,:] = allligandcoods[index,:]
	ocoods = np.float32(ocoods)

	tempdist = MDAnalysis.lib.distances.distance_array(ocoods,allligandcoods)
	A = np.argsort(tempdist)
	temp = int(A[0,1])

	Omatecood = np.zeros((1,3), dtype = float)
	Omatecood[0,:] = allligandcoods[temp,:]
	Omatecood = np.float32(Omatecood)

	tempdist2 = MDAnalysis.lib.distances.distance_array(Omatecood, allligandcoods)
	B = np.argsort(tempdist2)
	B = np.delete(B,0,axis = 1)
	for i in xrange(0,B.size):
		if B[0,i] == index:
			C = np.delete(B,i,axis = 1)
			break

	base1 = int(C[0,0])
	base2 = int(C[0,1])
	type1 = allligand[base1].type
	type2 = allligand[base2].type

	if type1 == 'O' or type2 == 'O':
		atype = 'carboxyl'
	else:
		atype = 'carbonyl'

	return atype

######################################################################################################################################################################
######################################################################################################################################################################

def Otypefinder(allligand,index,bond_dist):

	allligandcoods = allligand.positions
	ocoods = np.zeros((1,3), dtype = float)
	ocoods[0,:] = allligandcoods[index,:]
	ocoods = np.float32(ocoods)

	otype = 'NONE'

	tempdist = MDAnalysis.lib.distances.distance_array(allligandcoods, ocoods)
	A = np.where((tempdist < bond_dist) & (tempdist > 0.1))
	mates = np.ravel_multi_index(A, tempdist.shape)
	nummates = np.size(mates)

	numS = 0
	numP = 0
	numC = 0
	numN = 0
	numH = 0

	for j in mates:

		if allligand[j].type == 'C':
			numC = numC + 1
		elif allligand[j].type == 'P':
			numP = numP + 1

		elif allligand[j].type == 'N':
			numN = numN + 1
		elif allligand[j].type == 'H':
			numH = numH + 1
		elif allligand[j].type == 'S':
			numS = numS + 1

	if nummates == 1:
		if numC == 1:
			otype = carbonylorcarboxyl(allligand, index, bond_dist)
		elif numS == 1:
			otype = 'sulfone'
		elif numP == 1:
			otype = 'phosphone'
		elif numN == 1:
			otype = 'nitro'

	elif nummates == 2:
		if numC == 2:
			otype = 'ether'
		elif numH == 0:
			otype = 'ether'
		elif numH == 1:
			otype = 'hydroxyl'
		elif numP == 1:
			otype = 'phosphone'
		elif numS == 1:
			otype = 'sulfone'
	
	return otype
######################################################################################################################################################################
######################################################################################################################################################################
def Ntypefinder(allligand,index,bond_dist):
	allligandcoods = allligand.positions
	ncoods = np.zeros((1,3), dtype = float)
	ncoods[0,:] = allligandcoods[index,:]
	ncoods = np.float32(ncoods)

	tempdist = MDAnalysis.lib.distances.distance_array(allligandcoods, ncoods)
	A = np.where((tempdist < bond_dist) & (tempdist > 0.1))
	mates = np.ravel_multi_index(A, tempdist.shape)
	nummates = np.size(mates)
	
	ntype = 'NONE'

	numC = 0
	numN = 0
	numH = 0

	for j in mates:

		if allligand[j].type == 'C':
			numC = numC + 1
		elif allligand[j].type == 'N':
			numN = numN + 1
		elif allligand[j].type == 'H':
			numH = numH + 1
	

	if nummates == 1:
		if numC == 1:
			ntype = 'nitrile'

	elif nummates == 2:
		if numH == 0:
			ntype = 'imine'

	elif nummates == 3:

		if numH == 0:
			ntype = '3tertamine'
		elif numH == 1:
			ntype = '3secamine'
		elif numH == 2:
			ntype = '3priamine'

	elif nummates == 4:
		if numH == 1:
			ntype = 'secamine'
		elif numH == 2:
			ntype = 'priamine'
		elif numH == 3:
			ntype = 'ammonia'
		elif numH == 0:
			ntype = '4tertamine'


	return ntype
######################################################################################################################################################################
######################################################################################################################################################################
def carbonylwaters(allligand,index,bond_dist):

	D = 3.0

	allligandcoods = allligand.positions
	ocoods = np.zeros((1,3), dtype = float)
	ocoods[0,:] = allligandcoods[index,:]
	ocoods = np.float32(ocoods)

	tempdist = MDAnalysis.lib.distances.distance_array(ocoods,allligandcoods)
	A = np.argsort(tempdist)
	temp = int(A[0,1])
	Omatecood = np.zeros((1,3), dtype = float)
	Omatecood[0,:] = allligandcoods[temp,:]
	Omatecood = np.float32(Omatecood)

	tempdist2 = MDAnalysis.lib.distances.distance_array(Omatecood, allligandcoods)
	B = np.argsort(tempdist2)
	B = np.delete(B,0,axis = 1)
	for i in xrange(0,B.size):
		if B[0,i] == index:
			C = np.delete(B,i,axis = 1)

	base1cood = allligandcoods[C[0,0],:].copy()
	base2cood = allligandcoods[C[0,1],:].copy()

	v1 = Omatecood - base1cood
	v2 = Omatecood - base2cood

	xaxis = unitvector(np.cross(v1,v2))
	yaxis = unitvector(ocoods - Omatecood)
	zaxis = unitvector(np.cross(xaxis,yaxis))
	watercood = np.zeros((11,3), dtype = float)

	V = (D*yaxis)

	for i in xrange(0,5):
		degreeangle = -70.0 + (i*35)
		radangle = np.pi*(degreeangle/180)
		#Rotation Matrix
		#Vrot = vcos(@) + (k x v)sin(@) + k(k.v)(1 - cos(@))
		#Euler-Rodrigues Rotation Formula 
		#v - Vector to be rotated
		#@ - Angle to be rotated by
		#k - Unit Vector to be rotated about

		partA = V * np.cos(radangle)
		partB = np.sin(radangle) * np.cross(zaxis,V)
		partC = (1-np.cos(radangle)) * np.dot(zaxis,np.transpose(V)) * zaxis
		watercood[i,:] = partA + partB + partC

	for i in xrange(5,8):
		degreeangle = -35.0
		radangle = np.pi*(degreeangle/180)
		#Rotation Matrix
		#Vrot = vcos(@) + (k x v)sin(@) + k(k.v)(1 - cos(@))
		#Euler-Rodrigues Rotation Formula 
		#v - Vector to be rotated
		#@ - Angle to be rotated by
		#k - Unit Vector to be rotated about

		partX = watercood[i-4,:] * np.cos(radangle)
		partY = np.sin(radangle) * np.cross(xaxis,V)
		partZ = (1-np.cos(radangle)) * np.dot(xaxis,np.transpose(watercood[i-4,:])) * xaxis
		watercood[i,:] = partX + partY + partZ

	for i in xrange(8,11):
		degreeangle = 35.0
		radangle = np.pi*(degreeangle/180)
		#Rotation Matrix
		#Vrot = vcos(@) + (k x v)sin(@) + k(k.v)(1 - cos(@))
		#Euler-Rodrigues Rotation Formula 
		#v - Vector to be rotated
		#@ - Angle to be rotated by
		#k - Unit Vector to be rotated about

		partX = watercood[i-7,:] * np.cos(radangle)
		partY = np.sin(radangle) * np.cross(xaxis,V)
		partZ = (1-np.cos(radangle)) * np.dot(xaxis,np.transpose(watercood[i-7,:])) * xaxis
		watercood[i,:] = partX + partY + partZ
		

	watercood = watercood + ocoods

	return watercood
######################################################################################################################################################################
######################################################################################################################################################################
def carboxylwaters(allligand,index,bond_dist):

	D = 3.0
	
	allligandcoods = allligand.positions
	ocoods = np.zeros((1,3), dtype = float)
	ocoods[0,:] = allligandcoods[index,:]
	ocoods = np.float32(ocoods)

	tempdist = MDAnalysis.lib.distances.distance_array(ocoods,allligandcoods)
	A = np.argsort(tempdist)
	temp = int(A[0,1])
	Omatecood = np.zeros((1,3), dtype = float)
	Omatecood[0,:] = allligandcoods[temp,:]
	Omatecood = np.float32(Omatecood)

	tempdist2 = MDAnalysis.lib.distances.distance_array(Omatecood, allligandcoods)
	B = np.argsort(tempdist2)
	B = np.delete(B,0,axis = 1)
	for i in xrange(0,B.size):
		if B[0,i] == index:
			C = np.delete(B,i,axis = 1)

	base1cood = allligandcoods[C[0,0],:].copy()
	base2cood = allligandcoods[C[0,1],:].copy()

	v1 = Omatecood - base1cood
	v2 = Omatecood - base2cood

	xaxis = unitvector(np.cross(v1,v2))
	yaxis = unitvector(ocoods - Omatecood)
	zaxis = unitvector(np.cross(xaxis,yaxis))
	watercood = np.zeros((25,3), dtype = float)

	V = (D*yaxis)
	count = 0

	for i in xrange(0,5):
		
		degreeangle1 = -70.0 + (i*35)
		radangle1 = np.pi*(degreeangle1/180)
		#Rotation Matrix
		#Vrot = vcos(@) + (k x v)sin(@) + k(k.v)(1 - cos(@))
		#Euler-Rodrigues Rotation Formula 
		#v - Vector to be rotated
		#@ - Angle to be rotated by
		#k - Unit Vector to be rotated about

		partA = V * np.cos(radangle1)
		partB = np.sin(radangle1) * np.cross(zaxis,V)
		partC = (1-np.cos(radangle1)) * np.dot(zaxis,np.transpose(V)) * zaxis
		V2 = partA + partB + partC
		
		
		for j in xrange(0,5):

			degreeangle2 = -70.0 + (j*35)
			radangle2 = np.pi*(degreeangle2/180)
			#Rotation Matrix
			#Vrot = vcos(@) + (k x v)sin(@) + k(k.v)(1 - cos(@))
			#Euler-Rodrigues Rotation Formula 
			#v - Vector to be rotated
			#@ - Angle to be rotated by
			#k - Unit Vector to be rotated about

			partP = V2 * np.cos(radangle2)
			partQ = np.sin(radangle2) * np.cross(xaxis,V2)
			partR = (1-np.cos(radangle2)) * np.dot(xaxis,np.transpose(V2)) * xaxis
			watercood[count,:] = partP + partQ + partR
			count = count + 1


	watercood = watercood + ocoods

	return watercood
######################################################################################################################################################################
######################################################################################################################################################################
def secaminewater(allligand,index,bond_dist):

	D = 3.0

	allligandcoods = allligand.positions
	ncoods = np.zeros((1,3), dtype = float)
	ncoods[0,:] = allligandcoods[index,:]
	ncoods = np.float32(ncoods)

	tempdist = MDAnalysis.lib.distances.distance_array(allligandcoods, ncoods)
	A = np.where((tempdist < bond_dist) & (tempdist > 0.1))
	mates = np.ravel_multi_index(A, tempdist.shape)
	nummates = np.size(mates)
	hcoods = np.zeros((1,3), dtype = float)
	q = 0
	for j in mates:
		if allligand[j].type == 'H':
			hcoods[0,:] = allligandcoods[j,:]
			break

	watercood = np.zeros((1,3), dtype = float)

	hcoods = np.float32(hcoods)
	vector = unitvector(hcoods - ncoods)
	watercood[0,:] = ncoods + (D * vector)

	return watercood
######################################################################################################################################################################
######################################################################################################################################################################
def priaminewater(allligand,index,bond_dist):

	D = 3.0

	allligandcoods = allligand.positions
	ncoods = np.zeros((1,3), dtype = float)
	ncoods[0,:] = allligandcoods[index,:]
	ncoods = np.float32(ncoods)

	tempdist = MDAnalysis.lib.distances.distance_array(allligandcoods, ncoods)
	A = np.where((tempdist < bond_dist) & (tempdist > 0.1))
	mates = np.ravel_multi_index(A, tempdist.shape)
	nummates = np.size(mates)
	hcoods = np.zeros((2,3), dtype = float)

	i = 0
	for j in mates:
		if allligand[j].type == 'H':
			hcoods[i,:] = allligandcoods[j,:]
			i = i + 1

	hcoods = np.float32(hcoods)
	tempvector = hcoods - ncoods
	vector1 = unitvector(tempvector[0,:])
	vector2 = unitvector(tempvector[1,:])
	watercood = np.zeros((2,3), dtype = float)
	watercood[0,:] = ncoods + (D*vector1)
	watercood[1,:] = ncoods + (D*vector2)
	return watercood
######################################################################################################################################################################
######################################################################################################################################################################
def threetertaminewater(allligand,index,bond_dist):

	D = 3.0

	allligandcoods = allligand.positions
	ncoods = np.zeros((1,3), dtype = float)
	ncoods[0,:] = allligandcoods[index,:]
	ncoods = np.float32(ncoods)

	tempdist = MDAnalysis.lib.distances.distance_array(allligandcoods,ncoods)
	A = np.argsort(tempdist[:,0])

	base1cood = allligandcoods[A[1],:].copy()
	base2cood = allligandcoods[A[2],:].copy()
	base3cood = allligandcoods[A[3],:].copy()

	centroid = (base1cood + base2cood + base3cood)/3
	vector = unitvector(ncoods - centroid)

	watercood = np.zeros((1,3), dtype = float)
	watercood[0,:] = ncoods + (D*vector)
	return watercood
######################################################################################################################################################################
######################################################################################################################################################################
def threesecaminewater(allligand,index,bond_dist):

	D = 3.0

	allligandcoods = allligand.positions
	ncoods = np.zeros((1,3), dtype = float)
	ncoods[0,:] = allligandcoods[index,:]
	ncoods = np.float32(ncoods)

	watercood = np.zeros((2,3), dtype = float)

	tempdist = MDAnalysis.lib.distances.distance_array(allligandcoods,ncoods)
	A = np.argsort(tempdist[:,0])

	for i in xrange(1,4):
		if allligand[A[i]].type == 'H':
			vector1 = unitvector(allligandcoods[A[i],:] - ncoods[0,:])
			watercood[0,:] = ncoods + (D*vector1)

	base1cood = allligandcoods[A[1],:].copy()
	base2cood = allligandcoods[A[2],:].copy()
	base3cood = allligandcoods[A[3],:].copy()

	centroid = (base1cood + base2cood + base3cood)/3
	vector2 = unitvector(ncoods - centroid)
	watercood[1,:] = ncoods + (D*vector2)
	
	return watercood
######################################################################################################################################################################
######################################################################################################################################################################
def threepriaminewater(allligand,index,bond_dist):

	D = 3.0

	allligandcoods = allligand.positions
	ncoods = np.zeros((1,3), dtype = float)
	ncoods[0,:] = allligandcoods[index,:]
	ncoods = np.float32(ncoods)

	watercood = np.zeros((3,3), dtype = float)
	tempdist = MDAnalysis.lib.distances.distance_array(allligandcoods,ncoods)
	A = np.argsort(tempdist[:,0])

	j = 0
	for i in xrange(1,4):
		if allligand[A[i]].type == 'H':
			vector1 = unitvector(allligandcoods[A[i],:] - ncoods[0,:])
			watercood[j,:] = ncoods + (D*vector1)
			j = j + 1

	base1cood = allligandcoods[A[1],:].copy()
	base2cood = allligandcoods[A[2],:].copy()
	base3cood = allligandcoods[A[3],:].copy()

	centroid = (base1cood + base2cood + base3cood)/3
	vector2 = unitvector(ncoods - centroid)
	watercood[2,:] = ncoods + (D*vector2)
	
	return watercood
######################################################################################################################################################################
######################################################################################################################################################################
def ammoniawater(allligand,index,bond_dist):

	D = 3.0

	allligandcoods = allligand.positions
	ncoods = np.zeros((1,3), dtype = float)
	ncoods[0,:] = allligandcoods[index,:]
	ncoods = np.float32(ncoods)

	tempdist = MDAnalysis.lib.distances.distance_array(allligandcoods, ncoods)
	A = np.where((tempdist < bond_dist) & (tempdist > 0.1))
	mates = np.ravel_multi_index(A, tempdist.shape)
	nummates = np.size(mates)
	hcoods = np.zeros((3,3), dtype = float)

	i = 0
	for j in mates:
		if allligand[j].type == 'H':
			hcoods[i,:] = allligandcoods[j,:]
			i = i + 1

	hcoods = np.float32(hcoods)
	tempvector = hcoods - ncoods
	vector1 = unitvector(tempvector[0,:])
	vector2 = unitvector(tempvector[1,:])
	vector3 = unitvector(tempvector[2,:])
	watercood = np.zeros((3,3), dtype = float)
	watercood[0,:] = ncoods + (D*vector1)
	watercood[1,:] = ncoods + (D*vector2)
	watercood[2,:] = ncoods + (D*vector3)


	return watercood
######################################################################################################################################################################
######################################################################################################################################################################
def largeatom(allligand,index,bond_dist):

	D = 3.0
	
	allligandcoods = allligand.positions
	ocoods = np.zeros((1,3), dtype = float)
	ocoods[0,:] = allligandcoods[index,:]
	ocoods = np.float32(ocoods)

	tempdist = MDAnalysis.lib.distances.distance_array(ocoods,allligandcoods)
	A = np.argsort(tempdist)
	temp = int(A[0,1])
	Omatecood = np.zeros((1,3), dtype = float)
	Omatecood[0,:] = allligandcoods[temp,:]
	Omatecood = np.float32(Omatecood)

	tempdist2 = MDAnalysis.lib.distances.distance_array(Omatecood, allligandcoods)
	B = np.argsort(tempdist2)
	B = np.delete(B,0,axis = 1)
	for i in xrange(0,B.size):
		if B[0,i] == index:
			C = np.delete(B,i,axis = 1)

	base1cood = allligandcoods[C[0,0],:].copy()
	base2cood = allligandcoods[C[0,1],:].copy()

	v1 = base1cood - Omatecood
	v2 = base2cood - Omatecood
	tempvec = (v1 + v2)/2.0
	yaxis = unitvector(ocoods - Omatecood)
	zaxis = unitvector(np.cross(tempvec, yaxis))

	xaxis = unitvector(np.cross(yaxis,zaxis))

	
	watercood = np.zeros((25,3), dtype = float)

	V = (D*yaxis)
	count = 0

	for i in xrange(0,5):
		
		degreeangle1 = -70.0 + (i*35)
		radangle1 = np.pi*(degreeangle1/180)
		#Rotation Matrix
		#Vrot = vcos(@) + (k x v)sin(@) + k(k.v)(1 - cos(@))
		#Euler-Rodrigues Rotation Formula 
		#v - Vector to be rotated
		#@ - Angle to be rotated by
		#k - Unit Vector to be rotated about

		partA = V * np.cos(radangle1)
		partB = np.sin(radangle1) * np.cross(zaxis,V)
		partC = (1-np.cos(radangle1)) * np.dot(zaxis,np.transpose(V)) * zaxis
		V2 = partA + partB + partC
		
		
		for j in xrange(0,5):

			degreeangle2 = -70.0 + (j*35)
			radangle2 = np.pi*(degreeangle2/180)
			#Rotation Matrix
			#Vrot = vcos(@) + (k x v)sin(@) + k(k.v)(1 - cos(@))
			#Euler-Rodrigues Rotation Formula 
			#v - Vector to be rotated
			#@ - Angle to be rotated by
			#k - Unit Vector to be rotated about

			partP = V2 * np.cos(radangle2)
			partQ = np.sin(radangle2) * np.cross(xaxis,V2)
			partR = (1-np.cos(radangle2)) * np.dot(xaxis,np.transpose(V2)) * xaxis
			watercood[count,:] = partP + partQ + partR

			count = count + 1


	watercood = watercood + ocoods

	return watercood
######################################################################################################################################################################
######################################################################################################################################################################
def etherwater(allligand,index,bond_dist):

	D = 3.0

	allligandcoods = allligand.positions
	ocoods = np.zeros((1,3), dtype = float)
	ocoods[0,:] = allligandcoods[index,:]
	ocoods = np.float32(ocoods)

	tempdist = MDAnalysis.lib.distances.distance_array(ocoods,allligandcoods)
	A = np.argsort(tempdist)
	base1cood = np.zeros((1,3), dtype = float)
	base2cood = np.zeros((1,3), dtype = float)

	base1cood[0,:] = allligandcoods[A[0,1],:]
	base2cood[0,:] = allligandcoods[A[0,2],:]
	v1 = ocoods - base1cood
	v2 = ocoods - base2cood

	xaxis = np.cross(v1,v2)
	centroid = (base1cood + base2cood)/2
	yaxis = unitvector(ocoods - centroid)
	zaxis = unitvector(np.cross(xaxis,yaxis))

	watercood = np.zeros((5,3), dtype = float)
	V = D * yaxis

	for i in xrange(0,5):
		degreeangle = -70.0 + (i*35)
		radangle = np.pi*(degreeangle/180)
		#Rotation Matrix
		#Vrot = vcos(@) + (k x v)sin(@) + k(k.v)(1 - cos(@))
		#Euler-Rodrigues Rotation Formula 
		#v - Vector to be rotated
		#@ - Angle to be rotated by
		#k - Unit Vector to be rotated about

		partA = V * np.cos(radangle)
		partB = np.sin(radangle) * np.cross(zaxis,V)
		partC = (1-np.cos(radangle)) * np.dot(zaxis,np.transpose(V)) * zaxis
		watercood[i,:] = partA + partB + partC

	watercood = watercood + ocoods

	return watercood
######################################################################################################################################################################
######################################################################################################################################################################
def iminewater(allligand,index,bond_dist):

	D = 3.0

	allligandcoods = allligand.positions
	ncoods = np.zeros((1,3), dtype = float)
	ncoods[0,:] = allligandcoods[index,:]
	ncoods = np.float32(ncoods)

	tempdist = MDAnalysis.lib.distances.distance_array(ncoods,allligandcoods)
	A = np.argsort(tempdist)
	base1cood = np.zeros((1,3), dtype = float)
	base2cood = np.zeros((1,3), dtype = float)

	base1cood[0,:] = allligandcoods[A[0,1],:]
	base2cood[0,:] = allligandcoods[A[0,2],:]
	v1 = ncoods - base1cood
	v2 = ncoods - base2cood

	xaxis = np.cross(v1,v2)
	centroid = (base1cood + base2cood)/2
	yaxis = unitvector(ncoods - centroid)
	zaxis = unitvector(np.cross(xaxis,yaxis))

	watercood = np.zeros((5,3), dtype = float)
	V = D * yaxis

	for i in xrange(0,5):
		degreeangle = -70.0 + (i*35)
		radangle = np.pi*(degreeangle/180)
		#Rotation Matrix
		#Vrot = vcos(@) + (k x v)sin(@) + k(k.v)(1 - cos(@))
		#Euler-Rodrigues Rotation Formula 
		#v - Vector to be rotated
		#@ - Angle to be rotated by
		#k - Unit Vector to be rotated about

		partA = V * np.cos(radangle)
		partB = np.sin(radangle) * np.cross(zaxis,V)
		partC = (1-np.cos(radangle)) * np.dot(zaxis,np.transpose(V)) * zaxis
		watercood[i,:] = partA + partB + partC

	watercood = watercood + ncoods

	return watercood
######################################################################################################################################################################
######################################################################################################################################################################
def hydroxylwater(allligand,index,bond_dist):

	D = 3.0

	allligandcoods = allligand.positions
	ocoods = np.zeros((1,3), dtype = float)
	ocoods[0,:] = allligandcoods[index,:]
	ocoods = np.float32(ocoods)

	tempdist = MDAnalysis.lib.distances.distance_array(ocoods,allligandcoods)
	A = np.argsort(tempdist)
	base1cood = np.zeros((1,3), dtype = float)
	base2cood = np.zeros((1,3), dtype = float)
	
	for i in xrange(0,3):
		if allligand[A[0,i]].type == 'H':
			hindex = A[0,i]

	base1cood[0,:] = allligandcoods[A[0,1],:]
	base2cood[0,:] = allligandcoods[A[0,2],:]
	v1 = ocoods - base1cood
	v2 = ocoods - base2cood

	zaxis = np.cross(v1,v2)
	centroid = (base1cood + base2cood)/2
	yaxis = unitvector(ocoods - centroid)
	xaxis = unitvector(np.cross(yaxis,zaxis))

	watercood = np.zeros((12,3), dtype = float)
	V = D * yaxis
	count = 0

	for i in xrange(0,5):
		degreeangle = -70.0 + (i*35)
		radangle = np.pi*(degreeangle/180)
		#Rotation Matrix
		#Vrot = vcos(@) + (k x v)sin(@) + k(k.v)(1 - cos(@))
		#Euler-Rodrigues Rotation Formula 
		#v - Vector to be rotated
		#@ - Angle to be rotated by
		#k - Unit Vector to be rotated about

		partA = V * np.cos(radangle)
		partB = np.sin(radangle) * np.cross(xaxis,V)
		partC = (1-np.cos(radangle)) * np.dot(xaxis,np.transpose(V)) * xaxis
		Vec1 = partA + partB + partC

		for j in xrange(0,2):

			degreeangle2 = -17.5 + (j*35.0)
			radangle2 = np.pi*(degreeangle2/180)
			partX = Vec1 * np.cos(radangle2)
			partY = np.sin(radangle2) * np.cross(zaxis,Vec1)
			partZ = (1-np.cos(radangle2)) * np.dot(zaxis,np.transpose(Vec1)) * zaxis
			watercood[count,:] = partX + partY + partZ
			count = count + 1

	
	watercood = watercood + ocoods
	watercood[10,0] = 0
	watercood[10,1] = 0
	watercood[10,2] = 0
	
	ohvector = unitvector(allligandcoods[hindex,:] - ocoods)
	watercood[10,:] = ocoods + (D*ohvector)

	temppos = V + ocoods
	watercood[11,:] = (temppos + watercood[10,:])/2.0
	
	return watercood
######################################################################################################################################################################
######################################################################################################################################################################
def halogenwater(allligand,index,bond_dist):

	D = 3.0

	allligandcoods = allligand.positions
	halcoods = np.zeros((1,3), dtype = float)
	halcoods[0,:] = allligandcoods[index,:]
	halcoods = np.float32(halcoods)

	tempdist = MDAnalysis.lib.distances.distance_array(halcoods,allligandcoods,)
	A = np.argsort(tempdist)

	matecood = np.zeros((1,3), dtype = float)
	matecood[0,:] = allligandcoods[A[0,1],:]
	matecood = np.float32(matecood)

	vector = unitvector(halcoods - matecood)
	watercood = halcoods + (D*vector)

	return watercood
######################################################################################################################################################################
######################################################################################################################################################################
#THE ACTUAL PART STARTS HERE :)

def main(ligandinputfilename):

	bond_dist = 1.7
	delangle = 40.0
	delangle = ((delangle*np.pi)/180)
	U = MDAnalysis.Universe(ligandinputfilename)
	allligand = U.select_atoms('all')
	heavyligand = U.select_atoms('not type H')
	allligandcoods = allligand.positions
	numatoms = allligandcoods.shape[0]
	heavyligandcoods = heavyligand.positions
	waters = np.zeros((1,3), dtype = float)
	j = 0

	f1 = open('waterdetails.txt','w')


	for i in xrange(0,numatoms):
		atype = 'NONE'
		atom = str(allligand[i].type)
		if atom == 'O':
			atype = Otypefinder(allligand,i,bond_dist)
	
		elif atom == 'N':
			atype = Ntypefinder(allligand,i,bond_dist)
	
		elif atom == 'Cl' or atom == 'F' or atom == 'Br':
			atype = 'halogen'
		
	
		if atype == 'carbonyl':
			tempwatercood = carbonylwaters(allligand, i, bond_dist)
			waters = np.concatenate((waters,tempwatercood), axis = 0)
			f1.write(str(i))
			f1.write('\t')
			f1.write('2')
			f1.write('\n')

		elif atype == 'carboxyl':
			tempwatercood = carboxylwaters(allligand, i, bond_dist)
			waters = np.concatenate((waters,tempwatercood), axis = 0)
			f1.write(str(i))
			f1.write('\t')
			f1.write('2')
			f1.write('\n')
			
		elif atype == 'nitro':
			tempwatercood = largeatom(allligand, i, bond_dist)
			waters = np.concatenate((waters,tempwatercood), axis = 0)
			f1.write(str(i))
			f1.write('\t')
			f1.write('3')
			f1.write('\n')

		elif atype == 'hydroxyl':
			tempwatercood = hydroxylwater(allligand, i, bond_dist)
			waters = np.concatenate((waters,tempwatercood), axis = 0)
			f1.write(str(i))
			f1.write('\t')
			f1.write('2')
			f1.write('\n')

		elif atype == 'ether':
			tempwatercood = etherwater(allligand, i, bond_dist)
			waters = np.concatenate((waters,tempwatercood), axis = 0)
			f1.write(str(i))
			f1.write('\t')
			f1.write('2')
			f1.write('\n')

		elif atype == 'sulfone':
			tempwatercood = largeatom(allligand, i, bond_dist)
			waters = np.concatenate((waters,tempwatercood), axis = 0)
			f1.write(str(i))
			f1.write('\t')
			f1.write('2')
			f1.write('\n')

		elif atype == 'phosphone':
			tempwatercood = largeatom(allligand, i, bond_dist)
			waters = np.concatenate((waters,tempwatercood), axis = 0)
			f1.write(str(i))
			f1.write('\t')
			f1.write('2')
			f1.write('\n')

		elif atype == 'nitrile':
			#tempwatercood = nitrilewater(allligand, i, bond_dist)
			#waters = np.concatenate((waters,tempwatercood), axis = 0)
			f1.write(str(i))
			f1.write('\t')
			f1.write('1')
			f1.write('\n')

		elif atype == 'imine':
			tempwatercood = iminewater(allligand, i, bond_dist)
			waters = np.concatenate((waters,tempwatercood), axis = 0)
			f1.write(str(i))
			f1.write('\t')
			f1.write('1')
			f1.write('\n')

		elif atype == 'secamine':
			tempwatercood = secaminewater(allligand, i, bond_dist)
			waters = np.concatenate((waters,tempwatercood), axis = 0)
			f1.write(str(i))
			f1.write('\t')
			f1.write('1')
			f1.write('\n')

		elif atype == 'priamine':
			tempwatercood = priaminewater(allligand, i, bond_dist)
			waters = np.concatenate((waters,tempwatercood), axis = 0)
			f1.write(str(i))
			f1.write('\t')
			f1.write('2')
			f1.write('\n')

		elif atype == 'ammonia':

			tempwatercood = ammoniawater(allligand, i, bond_dist)
			waters = np.concatenate((waters,tempwatercood), axis = 0)
			f1.write(str(i))
			f1.write('\t')
			f1.write('3')
			f1.write('\n')

		elif atype == '3tertamine':
			tempwatercood = threetertaminewater(allligand, i, bond_dist)
			waters = np.concatenate((waters,tempwatercood), axis = 0)
			f1.write(str(i))
			f1.write('\t')
			f1.write('1')
			f1.write('\n')
			
		elif atype == '3secamine':
			tempwatercood = threesecaminewater(allligand, i , bond_dist)
			waters = np.concatenate((waters,tempwatercood), axis = 0)
			f1.write(str(i))
			f1.write('\t')
			f1.write('1')
			f1.write('\n')

		elif atype == '3priamine':
			tempwatercood = threepriaminewater(allligand, i, bond_dist)
			waters = np.concatenate((waters,tempwatercood), axis = 0)
			f1.write(str(i))
			f1.write('\t')
			f1.write('2')
			f1.write('\n')

		elif atype == '4tertamine':
			f1.write(str(i))
			f1.write('\t')
			f1.write('0')
			f1.write('\n')

		elif atype == 'halogen':
			tempwatercood = halogenwater(allligand, i, bond_dist)
			waters = np.concatenate((waters,tempwatercood), axis = 0)
			f1.write(str(i))
			f1.write('\t')
			f1.write('1')
			f1.write('\n')

	waters = np.delete(waters,0, axis = 0)
	waters = np.float32(waters)

	f1.close()
	tempwatligdist = MDAnalysis.lib.distances.distance_array(waters,heavyligandcoods)
	watligdist = np.amin(tempwatligdist, axis = 1)
	indices1 = np.where(watligdist < 2.10)
	waters = np.delete(waters, indices1, axis = 0)

	finalwaters = waters.copy()

	writewaterfile('placedwaters.pdb', finalwaters)

######################################################################################################################################################################

if __name__ == '__main__':

	main(sys.argv[1])
