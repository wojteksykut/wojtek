from math import sqrt
import numpy as np
from numpy.linalg import svd
import os


trash_atoms = ["H", "F", "Br", "Cl", "I"]

choice = input("Do you want to process one file or more?. Answer one or more.\n")
if choice == 'more':
	dir_path = input("Please enter a path to your directory with .xyz files:\n")
	print("New directory with output files will be saved as:\t" + dir_path + '_updated\n')
	paths = os.listdir(dir_path)
else:
	path = input("Please enter a path to your file with extension .xyz:\n")
	print("Output file will be saved as:\t" + path[:-4] + '_updated.xyz\n')
option = input("Do you want to specify ring type? Answer yes or no.\n")
if option == 'yes':
	option2 = input("Type in the atoms that make up the ring that you are intrested in with a space as separator. Eg. C C C C C C\n")
	user_cycle = option2.split()
option3 = input("Do you want two additional poits 1 angstrem far from centre of each cycle? Answer yes or no.\n")


def distance(a1, a2):
	return sqrt( (a1.x-a2.x)**2 + (a1.y-a2.y)**2 + (a1.z-a2.z)**2 )

def cycle_centre(list_of_atoms):
	x = sum([atom.x for atom in list_of_atoms])/len(list_of_atoms)
	y = sum([atom.y for atom in list_of_atoms])/len(list_of_atoms)
	z = sum([atom.z for atom in list_of_atoms])/len(list_of_atoms)
	return [round(x, 6),round(y, 6),round(z, 6)]

def planeFit(points):
    """
    p, n = planeFit(points)

    Given an array, points, of shape (d,...)
    representing points in d-dimensional space,
    fit an d-dimensional plane to the points.
    Return a point, p, on the plane (the point-cloud centroid),
    and the normal, n.
    """

    points = np.reshape(points, (np.shape(points)[0], -1)) # Collapse trialing dimensions
    assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1], points.shape[0])
    ctr = points.mean(axis=1)
    x = points - ctr[:,np.newaxis]
    M = np.dot(x, x.T) # Could also use np.cov(x) here.
    return ctr, svd(M)[0][:,-1]

def compare_cycles(user_cycle, cycle):
    n = len(cycle)
    for c in range(n):
        template = [user_cycle[(j+c)%n] for j in range(n)]
        i = 0
        while cycle[i].name == template[i] or cycle[i].name[0] == template[i]:
            if i == n-1:
                return cycle
            i+=1

class Atom:

	def __init__(self, nr, list_of_data):
		self.nr = nr
		self.name = list_of_data[0]
		self.x = float(list_of_data[1])
		self.y = float(list_of_data[2])
		self.z = float(list_of_data[3])
		self.ring = None
		self.neighbours = []

	def add_neighbour(self, atom):
		self.neighbours.append(atom)

	def show(self):
		nei = []
		for neighbour in self.neighbours:
			nei.append((neighbour.nr,neighbour.name))
		print(self.nr, self.name, self.x, self.y, self.z, nei)


class Graph:

	def __init__(self, atoms, bonds = None):
		
		self.atoms = atoms
		self.bonds = bonds
		if self.bonds == None:
			self.bonds = []
			for count, a1 in enumerate(self.atoms):
				if a1.name[0] not in trash_atoms:
					for a2 in self.atoms[count+1:]:
						if a2.name[0] not in trash_atoms:
							if distance(a1, a2) < 2:
								self.bonds.append((a1,a2))
								a1.add_neighbour(a2)
								a2.add_neighbour(a1)


	def find_cycles(self, root = None):
		
		gnodes = set(self.atoms)
		cycles = []
		while gnodes:
			if root == None:
				root = gnodes.pop()
			stack = [root]
			pred = {root: root}
			used = {root: set()}
			while stack:

				current = stack.pop()
				current_used = used[current]

				for neighbour in current.neighbours:

					if neighbour not in used:
						pred[neighbour] = current
						stack.append(neighbour)
						used[neighbour] = {current}

					elif neighbour not in current_used:
						pn = used[neighbour]
						cycle = [neighbour, current]

						p = pred[current]

						while p not in pn:
							cycle.append(p)
							p = pred[p]
						cycle.append(p)
						cycles.append(cycle)
						used[neighbour].add(current)
			gnodes -= set(pred)
			root = None
		return cycles



if choice == 'one':

	file = open(path, "r").readlines()
	atoms = []
	for count, line in enumerate(file[2:],1):
		if line.split() != []:
			atoms.append(Atom(count, line.split()))

	G = Graph(atoms)
	C = G.find_cycles()

	if option == 'yes':
		C_user = []
		for cycle in C:
			print('cykl: \t', [atom.name for atom in cycle])
			if len(cycle) == len(user_cycle):
				com_cycle = compare_cycles(user_cycle, cycle)
				if com_cycle != None:
					C_user.append(com_cycle)
					#print('C_user: \t', [[atom.name for atom in cycle] for cycle in C_user])

	if option == 'yes':
		C = C_user
	num_of_cycles = len(C)

	input("There are %d cycles found. Press enter to close" % (num_of_cycles))

	new_path = path[:-4] + '_updated.xyz'
	file2 = open(new_path, "w")

	for count, line in enumerate(file):
	    if count == 0:
	    	if option3 == 'yes':
	    		file2.write(str(int(line.split()[0]) + num_of_cycles*3)+'\n')
	    	else:
	        	file2.write(str(int(line.split()[0]) + num_of_cycles)+'\n')
	    else:
	        if line.split() != []:
	            file2.write(line)

	if option == 'yes':
		C = C_user

	if C != [None]:
		for cycle in C:
			S = cycle_centre(cycle)
			file2.write('Bq' + '\t' + str(S[0]) + '\t'+ str(S[1]) + '\t'+ str(S[2]) + '\n')
			
			if option3 == 'yes':
				X,Y,Z = [],[],[]
				for atom in cycle:
					X.append(atom.x)
					Y.append(atom.y)
					Z.append(atom.z)
				N = planeFit([X,Y,Z])[1]

				file2.write('Bq' + '\t' + str(round(S[0]+N[0],6)) + '\t'+ str(round(S[1]+N[1],6)) + '\t'+ str(round(S[2]+N[2],6)) + '\n')
				file2.write('Bq' + '\t' + str(round(S[0]-N[0],6)) + '\t'+ str(round(S[1]-N[1],6)) + '\t'+ str(round(S[2]-N[2],6)) + '\n')
	file2.close()

if choice == 'more':
	
	new_dir = dir_path + '_updated'
	os.mkdir(new_dir)
	for path in paths:
		file = open(dir_path + '\\' + path, "r").readlines()
		atoms = []
		for count, line in enumerate(file[2:],1):
			if line.split() != []:
				atoms.append(Atom(count, line.split()))

		G = Graph(atoms)
		C = G.find_cycles()

		if option == 'yes':
			C_user = []
			for cycle in C:
				if len(cycle) == len(user_cycle):
					com_cycle = compare_cycles(user_cycle, cycle)
					if com_cycle != None:
						C_user.append(com_cycle)
						
		if option == 'yes':
			C = C_user
		num_of_cycles = len(C)

		new_path = new_dir + '\\' + path[:-4] + '_updated.xyz'
		file2 = open(new_path, "w")

		for count, line in enumerate(file):
		    if count == 0:
		    	if option3 == 'yes':
		    		file2.write(str(int(line.split()[0]) + num_of_cycles*3)+'\n')
		    	else:
		        	file2.write(str(int(line.split()[0]) + num_of_cycles)+'\n')
		    else:
		        if line.split() != []:
		            file2.write(line)

		if option == 'yes':
			C = C_user

		if C != [None]:
			for cycle in C:
				S = cycle_centre(cycle)
				file2.write('Bq' + '\t' + str(S[0]) + '\t'+ str(S[1]) + '\t'+ str(S[2]) + '\n')
				
				if option3 == 'yes':
					X,Y,Z = [],[],[]
					for atom in cycle:
						X.append(atom.x)
						Y.append(atom.y)
						Z.append(atom.z)
					N = planeFit([X,Y,Z])[1]
					file2.write('Bq' + '\t' + str(round(S[0]+N[0],6)) + '\t'+ str(round(S[1]+N[1],6)) + '\t'+ str(round(S[2]+N[2],6)) + '\n')
					file2.write('Bq' + '\t' + str(round(S[0]-N[0],6)) + '\t'+ str(round(S[1]-N[1],6)) + '\t'+ str(round(S[2]-N[2],6)) + '\n')
		file2.close()


