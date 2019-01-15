####---------------imports------------


import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter


####---------------functions---------------


def parsing_file(filename):
	'''
	returns a list with all the lines of the file
	'''
	print("Parsing file "+filename) 
	f=open(filename)
	f_lines=f.readlines() ##save the entire file into an array
	f.close()
	return f_lines

def file_matrices(a_list):
	genotypes=[]
	header = a_list[0].split(" ")
	nsam=int(header[1])
	nreps=int(header[2])
	x_ton_positions = {i:[] for i in range(1,nsam)}
	for i in range (5,len(a_list),nsam+4):
		
		###positions
		
		genotypes_pos=np.array(a_list[i].split(" ")[1:-1], dtype=float)
		pos_left=genotypes_pos[genotypes_pos<=0.5]
		pos_right=genotypes_pos[genotypes_pos>0.5]

		
		genotypes=np.array([list(a_list[i+q])[0:-1] for q in range(1,nsam+1)], dtype=int)
		left=np.array(genotypes, dtype=int)[:,genotypes_pos<=0.5].sum(axis=0)
		right=np.array(genotypes, dtype=int)[:,genotypes_pos>0.5].sum(axis=0)
		
		for q in range(len(pos_left)):
			x_ton_positions[left[q]].append(pos_left[q])
			
		for q in range(len(pos_right)):
			x_ton_positions[right[q]].append(1-pos_right[q])
	return x_ton_positions		
	
		
		
####---------------main----------------------


for f in range(1,61):
	
	file_list=parsing_file("bottleneck_simulation_results_only/ms_neutral{}.out".format(str(f)))
	
	selection=file_matrices(file_list)
	# An "interface" to matplotlib.axes.Axes.hist() method
	n, bins, patches = plt.hist(x=[selection[1],selection[2],selection[18],selection[19]], color=['blue','green','red','yellow'], alpha=0.7, rwidth=0.85)
	plt.grid(axis='y', alpha=0.75)
	plt.xlabel('Distance from selection')
	plt.ylabel('p(x=K)')
	plt.legend(["K=1","K=2","K=18",'K=19'])
	plt.title('1D SFS Histogram')

	# Set a clean upper y-axis limit.

	plt.savefig('1_d_distance_plots/neutral_{}.png'.format(str(f)))
	plt.clf()
	plt.cla()
	plt.close()

		

	