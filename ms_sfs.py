###ask if the ms file needs cleaning

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
############------run ms from python------


###-----------------functions----------------

def ms_simulation(rt,s, t , r):
	#os.system("ms "+n+" "+k+" -t "+t+" -r "+r+" 500 > ms_"+t+"_"+r+"_500 -eG 0.5 2") ####-eN 0.1
	#os.system("grep -v '[^01]' ms_"+t+"_"+r+"_500  > ms_"+t+"_"+r+".clean")
	return None
	
def compare_matrices(matrix1, matrix2) :
	
	assert(matrix1.shape==matrix2)
	q= matrix1 -matrix2
	print (q)
	

def sfs_matrix(vector, matrix) :
	'''
	returns a 2d matrix with rows as an SNP count and columns as the immediate SNP count of the next site 
	'''
	
	##print (len(vector))
	for i in range(0,len(vector)-1) :
		matrix[int(vector[i]-1),int(vector[i+1]-1)]+=1
	return matrix/len(vector)

def sfs_vector(n, vector):
	sfs=[]
	for i in range(1,n) :
		sfs.append(np.count_nonzero(vector == i))
	return sfs	

def matrix_cluster(matrix) :
	for i in range(matrix.shape[0]) :
		for q in range(matrix.shape[1]) :
			pass
		
##def open_file(filename):
##	f=open(filename)
##	f_lines=f.readlines() ##save the entire file into an array
##	f.close()
##	return f_lines
	
def matrix_heatmap(a_matrix, output):
	plt.imshow(a_matrix, cmap='viridis', interpolation='nearest')
	plt.savefig(output)
	return None
	
	
def genotypes_distance(a_list, dist, frame):
	header = a_list[0].split(" ")
	n = int(header[1])
	k = header[2]
	score_matrices=[]
	for i in range (5,len(lines),n+4) :  ##we know the number of samples since we gave it as an argument in the ms
		pos_header=np.array(lines[i].split(" ")[1:-1], dtype=float) 
####logiika tha mpei for gia to frame
		left=pos_header[pos_header<dist]
		right= pos_header[pos_header>1-dist]

		score_matrix = np.zeros((n-1,n-1)) ###maximum shape (n-1,n-1)
		number_of_snps = len(left)+len(right)
		genotype_array=np.zeros((n,number_of_snps))
		for g in range (0,n) :
			left_genos=list(lines[i+g+1])[0:len(left)]
			right_genos=list(lines[i+g+1])[len(pos_header)-len(right):-1]
			genotype_array[g] = left_genos + right_genos

		sth_tons=genotype_array.sum(axis=0)
		##sfs= sfs_vector(n, sth_tons)
		score_matrix=sfs_matrix(sth_tons, score_matrix)
		##print (sum(np.sum(score_matrix,axis=1)))
		##np.savetxt("test.txt",score_matrix.astype(int), fmt="%d")
	
		#assert(0)
		
		score_matrices.append(score_matrix)
	return score_matrices
##################-------main-------------

for fc in range(16,61) :	
	print("Parsing file number"+str(fc))
	lines=open_file("bottleneck_simulation_results_only/mssel{}.out".format(str(fc)))
	score_matrices1=genotypes_distance(lines, 0.2, 10)

	b=sum(score_matrices1)/len(score_matrices1)
	matrix_heatmap(b,"ms_selection_boxplot_{}_0.2.png".format(str(fc)))
	#assert (0)
#np.savetxt("file.txt",
#matrix_list=matrix_from_file("10","20", "100", "200" )
#)
#for i in matrix_list:cd

	#print (i)
	#plt.hist(i)
	#fig = plt.figure()
	#ax = fig.add_subplot(111, projection='3d')
	#xs = np.arange(10)
	#ys = np.arange(10)
	#ax.bar(xs, ys, zs=z, zdir='y', alpha=0.8)

	#plt.show()
	
####heatmap
####check zero clusters

		
