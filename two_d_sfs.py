####---------------imports------------


import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


####---------------functions---------------

	
def matrix_heatmap(a_matrix, output):
	plt.imshow(a_matrix, cmap='viridis', interpolation='nearest')
	plt.savefig(output)
	return None
	

def sfs_matrix(vector, n) :
	'''
	returns a 2d matrix with rows as an SNP count and columns as the immediate SNP count of the next site 
	'''
	##print (len(vector))
	matrix = np.zeros((n-1,n-1))
	for i in range(0,len(vector)-1) :
		matrix[int(vector[i]-1),int(vector[i+1]-1)]+=1
	return matrix/len(vector)


def parsing_file(filename):
	'''
	returns a list with all the lines of the file
	'''
	print("Parsing file "+filename)
	f=open(filename)
	f_lines=f.readlines() ##save the entire file into an array
	f.close()
	return f_lines

def list_to_arrays(a_list, start, window):
	genotypes=[]
	header = a_list[0].split(" ")
	nsam=int(header[1])
	nreps=int(header[2])
	total_score_matrix = np.zeros((nsam-1,nsam-1))
	for i in range (5,len(a_list),nsam+4):
		genotypes_pos=np.array(a_list[i].split(" ")[1:-1], dtype=float)

		genotypes=[list(a_list[i+q])[0:-1]for q in range(1,nsam+1)]
		left=np.array(genotypes, dtype=int)[:,np.bitwise_or(genotypes_pos>start,genotypes_pos<start+window)]
		right=np.array(genotypes, dtype=int)[:,np.bitwise_or(genotypes_pos<1-start,genotypes_pos>1-(start+window))]
		score_matrix=sfs_matrix(np.hstack((left,right)).sum(axis=0), nsam)
		total_score_matrix += score_matrix
	print(total_score_matrix/nreps)
	matrix_heatmap(total_score_matrix/nreps, "test.pdf")
	assert(0)
	

#		genotypes.append(genotypes_pos.astype(np.float))
#	print(genotypes)
####---------------main----------------------

file_list=parsing_file("bottleneck_simulation_results_only/mssel{}.out".format(str(60)))
list_to_arrays(file_list,0.2,0.2)

