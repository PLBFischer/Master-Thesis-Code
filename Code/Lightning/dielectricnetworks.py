import numpy as np

def create_round_domain(L):

	# Create a square lattice covering the entire circle
	N = 2*L + 1

	n = N*N
	m = (n)*2 - N - N

	# Generate the node-edge incidence matrix
	bondList = np.zeros([m,2],dtype=int) 
	nodeList = np.zeros([n,2],dtype=int) 
	nodeID = np.arange(0,n)

	for k in range(N):
		for j in range(N):
			idx = k*N + j
			nodeList[idx,0] = j
			nodeList[idx,1] = k            

	# to right
	l = 0
	for k in range(N):
		for j in range(N-1):
			idx = k*N + j
			bondList[l,0] =  idx
			bondList[l,1] =  idx+1                  
			l += 1

	# to top
	for k in range(N-1):
		for j in range(N):
			idx = k*N + j
			bondList[l,0] =  idx
			bondList[l,1] =  idx+N                      
			l += 1

	# Determine the distance of each node to the center
	R = nodeList.astype(float)
	R = R - np.mean(R,axis=0)

	dR = np.sqrt(np.sum(np.power(R,2),axis=1))

	# Remove nodes and bonds that are not within the domain of interesrt
	nodemask = dR < L
	dR = dR[nodemask]

	nodeList = nodeList[nodemask]
	nodeID = nodeID[nodemask]

	bondmask = np.isin(bondList, nodeID)
	bondmask = np.logical_and(bondmask[:,0],bondmask[:,1])
	bondList = bondList[bondmask]

	# Renumber nodes in bondList
	newnodeID = np.arange(0,len(nodeID),dtype=int)
	nodeID_mask = np.zeros(n,dtype=int)
	nodeID_mask[nodemask] = newnodeID
	bondList = nodeID_mask[bondList]

	# Determine source nodes
	sourcemask = np.isclose(dR,0.0)

	# determine boundary nodes
	_, counts = np.unique(bondList,return_counts=True)
	sinkmask = counts < 4

	# Remove any edges between boundary nodes
	boundarybondmask = np.isin(bondList, newnodeID[sinkmask])
	boundarybondmask = np.invert(np.logical_and(boundarybondmask[:,0],boundarybondmask[:,1]))	
	bondList = bondList[boundarybondmask]

	return bondList, nodeList, sourcemask, sinkmask

def create_pointsource_pointsink_circulardomain(L,d):


	# Create a square lattice covering the entire circle
	N = 2*L + 1

	n = N*N
	m = (n)*2 - N - N

	# Generate the node-edge incidence matrix
	bondList = np.zeros([m,2],dtype=int) 
	nodeList = np.zeros([n,2],dtype=int) 
	nodeID = np.arange(0,n)

	for k in range(N):
		for j in range(N):
			idx = k*N + j
			nodeList[idx,0] = j
			nodeList[idx,1] = k            

	# to right
	l = 0
	for k in range(N):
		for j in range(N-1):
			idx = k*N + j
			bondList[l,0] =  idx
			bondList[l,1] =  idx+1                  
			l += 1

	# to top
	for k in range(N-1):
		for j in range(N):
			idx = k*N + j
			bondList[l,0] =  idx
			bondList[l,1] =  idx+N                      
			l += 1

	# Determine the distance of each node to the center
	R = nodeList.astype(float)
	R = R - np.mean(R,axis=0)

	dR = np.sqrt(np.sum(np.power(R,2),axis=1))

	# Remove nodes and bonds that are not within the domain of interesrt
	nodemask = dR < L
	dR = dR[nodemask]
	R = R[nodemask]

	nodeList = nodeList[nodemask]
	nodeID = nodeID[nodemask]

	bondmask = np.isin(bondList, nodeID)
	bondmask = np.logical_and(bondmask[:,0],bondmask[:,1])
	bondList = bondList[bondmask]

	# Renumber nodes in bondList
	newnodeID = np.arange(0,len(nodeID),dtype=int)
	nodeID_mask = np.zeros(n,dtype=int)
	nodeID_mask[nodemask] = newnodeID
	bondList = nodeID_mask[bondList]

	# Determine source nodes
	sourcemask = np.logical_and(np.isclose(R[:,0],-0.5*d),np.isclose(R[:,1],0.0))

	# Determine sink nodes
	sinkmask = np.logical_and(np.isclose(R[:,0],0.5*d),np.isclose(R[:,1],0.0))

	return bondList, nodeList, sourcemask, sinkmask

def squaredomain(N):

	n = N*N
	m = (n)*2 - N - N

	# Generate the node-edge incidence matrix
	bondList = np.zeros([m,2],dtype=int) 
	nodeList = np.zeros([n,2],dtype=int) 
	nodeID = np.arange(0,n)

	for k in range(N):
		for j in range(N):
			idx = k*N + j
			nodeList[idx,0] = j
			nodeList[idx,1] = k            

	# to right
	l = 0
	for k in range(N):
		for j in range(N-1):
			idx = k*N + j
			bondList[l,0] =  idx
			bondList[l,1] =  idx+1                  
			l += 1

	# to top
	for k in range(N-1):
		for j in range(N):
			idx = k*N + j
			bondList[l,0] =  idx
			bondList[l,1] =  idx+N                      
			l += 1

	return bondList, nodeList, nodeID

def rectangulardomain(W, L):
    n = W * L
    m = (W - 1) * L + W * (L - 1)

    # Generate the node-edge incidence matrix
    bondList = np.zeros([m, 2], dtype=int)
    nodeList = np.zeros([n, 2], dtype=int)
    nodeID = np.arange(0, n)

    for k in range(L):
        for j in range(W):
            idx = k * W + j
            nodeList[idx, 0] = j
            nodeList[idx, 1] = k

    # to right
    l = 0
    for k in range(L):
        for j in range(W - 1):
            idx = k * W + j
            bondList[l, 0] = idx
            bondList[l, 1] = idx + 1
            l += 1

    # to top
    for k in range(L - 1):
        for j in range(W):
            idx = k * W + j
            bondList[l, 0] = idx
            bondList[l, 1] = idx + W
            l += 1

    return bondList, nodeList, nodeID

def create_pointsource_pointsink_squaredomain(L):

	# Create a square lattice covering the entire circle
	N = 2*L + 1

	bondList, nodeList, nodeID = squaredomain(N)

	# Determine the distance of each node to the center
	R = nodeList.astype(float)
	R = R - np.mean(R,axis=0)
 
	print("#####",R)


	# Determine source nodes
	sourcemask = np.logical_and(np.isclose(R[:,0],0),np.isclose(R[:,1],-L))

	# Determine sink nodes
	sinkmask = np.logical_and(np.isclose(R[:,0],0),np.isclose(R[:,1],L))

	return bondList, nodeList, sourcemask, sinkmask


def create_linesource_linesink_squaredomain(L):

	# Create a square lattice covering the entire circle
	N = 2*L + 1

	bondList, nodeList , nodeID = squaredomain(N)

	# Determine the distance of each node to the center
	R = nodeList.astype(float)
	R = R - np.mean(R,axis=0)

	# Determine source nodes
	sourcemask = np.isclose(R[:,1],-L)

	# Determine sink nodes
	sinkmask = np.isclose(R[:,1],L)

	# Remove any edges between sources nodes
	boundarybondmask = np.isin(bondList, nodeID[sourcemask])
	boundarybondmask = np.invert(np.logical_and(boundarybondmask[:,0],boundarybondmask[:,1]))	
	bondList = bondList[boundarybondmask]

	# Remove any edges between sink nodes
	boundarybondmask = np.isin(bondList, nodeID[sinkmask])
	boundarybondmask = np.invert(np.logical_and(boundarybondmask[:,0],boundarybondmask[:,1]))	
	bondList = bondList[boundarybondmask]

	return bondList, nodeList, sourcemask, sinkmask


def create_linesource_pointsink_squaredomain(L):

	# Create a square lattice covering the entire circle
	N = 2*L + 1

	bondList, nodeList , nodeID = squaredomain(N)

	# Determine the distance of each node to the center
	R = nodeList.astype(float)
	R = R - np.mean(R,axis=0)

	# Determine source nodes
	sourcemask = np.isclose(R[:,1],-L)

	# Determine sink nodes
	sinkmask = np.logical_and(np.isclose(R[:,0],0),np.isclose(R[:,1],L))

	# Remove any edges between sources nodes
	boundarybondmask = np.isin(bondList, nodeID[sinkmask])
	boundarybondmask = np.invert(np.logical_and(boundarybondmask[:,0],boundarybondmask[:,1]))	
	bondList = bondList[boundarybondmask]

	return bondList, nodeList, sourcemask, sinkmask


def create_pointsource_linesink_squaredomain(L):

	# Create a square lattice covering the entire circle
	N = 2*L + 1

	bondList, nodeList , nodeID = squaredomain(N)

	# Determine the distance of each node to the center
	R = nodeList.astype(float)
	R = R - np.mean(R,axis=0)

	# Determine source nodes
	sourcemask = np.logical_and(np.isclose(R[:,0],0),np.isclose(R[:,1],-L))

	# Determine sink nodes
	sinkmask = np.isclose(R[:,1],L)

	# Remove any edges between sources nodes
	boundarybondmask = np.isin(bondList, nodeID[sourcemask])
	boundarybondmask = np.invert(np.logical_and(boundarybondmask[:,0],boundarybondmask[:,1]))	
	bondList = bondList[boundarybondmask]

	return bondList, nodeList, sourcemask, sinkmask

def create_pointsource_linesink_rectangulardomain(W, L):

	NW = 2*W + 1
	NL = 2*L + 1

	bondList, nodeList , nodeID = rectangulardomain(NW, NL)

	# Determine the distance of each node to the center
	R = nodeList.astype(float)
	R = R - np.mean(R,axis=0)

	# Determine source nodes
	sourcemask = np.logical_and(np.isclose(R[:,0],0),np.isclose(R[:,1],-L))

	# Determine sink nodes
	sinkmask = np.isclose(R[:,1],L)

	# Remove any edges between sources nodes
	boundarybondmask = np.isin(bondList, nodeID[sourcemask])
	boundarybondmask = np.invert(np.logical_and(boundarybondmask[:,0],boundarybondmask[:,1]))	
	bondList = bondList[boundarybondmask]

	return bondList, nodeList, sourcemask, sinkmask



def create_pointsource_linesidesink_squaredomain(L):

	# Create a square lattice covering the entire circle
	N = 2*L + 1

	bondList, nodeList , nodeID = squaredomain(N)

	# Determine the distance of each node to the center
	R = nodeList.astype(float)
	R = R - np.mean(R,axis=0)

	# Determine source nodes
	sourcemask = np.logical_and(np.isclose(R[:,0],0),np.isclose(R[:,1],-L))

	# Determine sink nodes
	sinkmask = np.logical_or(np.isclose(R[:,1],L),np.isclose(np.abs(R[:,0]),L))

	# Remove any edges between sources nodes
	boundarybondmask = np.isin(bondList, nodeID[sourcemask])
	boundarybondmask = np.invert(np.logical_and(boundarybondmask[:,0],boundarybondmask[:,1]))	
	bondList = bondList[boundarybondmask]

	return bondList, nodeList, sourcemask, sinkmask

def create_pointsource_linepartsidesink_squaredomain(L,p):

	# Create a square lattice covering the entire circle
	N = 2*L + 1

	bondList, nodeList , nodeID = squaredomain(N)

	# Determine the distance of each node to the center
	R = nodeList.astype(float)
	R = R - np.mean(R,axis=0)

	# Determine source nodes
	sourcemask = np.logical_and(np.isclose(R[:,0],0),np.isclose(R[:,1],-L))

	# Determine sink nodes
	sinkmask = np.logical_or(np.isclose(R[:,1],L),np.logical_and(np.isclose(np.abs(R[:,0]),L),R[:,1]>p))

	# Remove any edges between sources nodes
	boundarybondmask = np.isin(bondList, nodeID[sourcemask])
	boundarybondmask = np.invert(np.logical_and(boundarybondmask[:,0],boundarybondmask[:,1]))	
	bondList = bondList[boundarybondmask]

	return bondList, nodeList, sourcemask, sinkmask

def triangular_lattice_point_point(rows,cols):
	bonds = []
	coors = []

	for i in range(rows):
		y0 = 1/4 * (1 + (-1)**i)
		for j in range(cols):
			x = i/2
			y = y0 + j
			coors.append(np.array([x,y]))

			index = cols * i + j
			dictn = (-1)**i

			if i != rows-1:
				bonds.append(np.array([index,index+cols]))
			if j != cols-1:
				bonds.append(np.array([index,index+1]))

			if j != 0 and dictn == -1 and i!=rows-1:
				bonds.append(np.array([index,index+cols-1]))
			if j != cols-1 and dictn == 1 and i!=rows-1:
				bonds.append(np.array([index,index+cols+1]))

	coors = 2*(np.array(coors))
	coors = coors.astype(int)
	bonds = np.array(bonds)

	source_mask = np.zeros(coors.shape[0])
	source_mask[int(cols/2)] = 1
	source_mask = source_mask.astype(bool)

	sink_mask = np.zeros(coors.shape[0])
	sink_mask[int((rows-1)*cols + cols/2)] = 1
	sink_mask = sink_mask.astype(bool)
	return bonds, coors, source_mask, sink_mask

def triangular_lattice_point_line(rows,cols):
	bonds = []
	coors = []

	for i in range(rows):
		y0 = 1/4 * (1 + (-1)**i)
		for j in range(cols):
			x = i/2
			y = y0 + j
			coors.append(np.array([x,y]))

			index = cols * i + j
			dictn = (-1)**i

			if i != rows-1:
				bonds.append(np.array([index,index+cols]))
			if j != cols-1:
				bonds.append(np.array([index,index+1]))

			if j != 0 and dictn == -1 and i!=rows-1:
				bonds.append(np.array([index,index+cols-1]))
			if j != cols-1 and dictn == 1 and i!=rows-1:
				bonds.append(np.array([index,index+cols+1]))

	coors = 2*(np.array(coors))
	coors = coors.astype(int)
	bonds = np.array(bonds)

	source_mask = np.zeros(coors.shape[0])
	source_mask[int(cols/2)] = 1
	source_mask = source_mask.astype(bool)

	sink_mask = np.zeros(coors.shape[0])
	sink_mask[(rows-1)*cols:] = 1
	sink_mask = sink_mask.astype(bool)
	return bonds, coors, source_mask, sink_mask


def triangular_lattice_line_line(rows,cols):
	bonds = []
	coors = []

	for i in range(rows):
		y0 = 1/4 * (1 + (-1)**i)
		for j in range(cols):
			x = i/2
			y = y0 + j
			coors.append(np.array([x,y]))

			index = cols * i + j
			dictn = (-1)**i

			if i != rows-1:
				bonds.append(np.array([index,index+cols]))
			if j != cols-1:
				bonds.append(np.array([index,index+1]))

			if j != 0 and dictn == -1 and i!=rows-1:
				bonds.append(np.array([index,index+cols-1]))
			if j != cols-1 and dictn == 1 and i!=rows-1:
				bonds.append(np.array([index,index+cols+1]))

	coors = 2*(np.array(coors))
	coors = coors.astype(int)
	bonds = np.array(bonds)

	source_mask = np.zeros(coors.shape[0])
	source_mask[:cols] = 1
	source_mask = source_mask.astype(bool)

	sink_mask = np.zeros(coors.shape[0])
	sink_mask[(rows-1)*cols:] = 1
	sink_mask = sink_mask.astype(bool)
	return bonds, coors, source_mask, sink_mask
