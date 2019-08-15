from networkx import grid_graph
import networkx as nx
import numpy as np
import random
from matplotlib import colors, cm
import matplotlib.pyplot as plt
import pandas as pd
#%matplotlib qt

import matplotlib.animation as animation
#from dendropy.calculate import treemeasure
#from dendropy import Tree as DTree

from random import sample, random, choice
#from copy import deepcopy
#from math import log
#from matplotlib. import PillowWriter
#from ete3 import Tree, NodeStyle, TreeStyle
#import tqdm
#from Bio import Phylo
#from plot_eteTree import plot_tree
from collections import Counter
#import pandas as pd
from scipy.stats import poisson 

import warnings
warnings.filterwarnings('ignore')

# Takes a genotype and converts it to an integer for use indexing the fitness landscape list 
def convertGenotypeToInt(genotype):
	out = 0
	for bit in genotype:
		out = (out << 1) | bit
	return out

# Converts an integer to a genotype by taking the binary value and padding to the left by 0s		
def convertIntToGenotype(anInt, pad):
	offset = 2**pad
	return [int(x) for x in bin(offset+anInt)[3:]]	

def flip(allele):
    if allele ==1:
        allele = 0
    else:
        allele =1
        
    return allele


def fitness(d,m,k,s_p,s,pk):
    p = s - s_p
    y = s + (p)* ((d**k)/((d**k)- (s-p)/s))*((m - pk)/pk)
    return y


def mutation(r,dist):
    check = np.append(dist,r)
    check.sort()
    m = np.where(check ==r )[0][0]
    return m
    
def poisson_max_cdf(x,mu,n):
    y = poisson.cdf(x,mu)**n
    return y
    
    
def sort_pairs(pair):
    # Extract integer after "r".
    return int(pair[0][1:])
   
def make_tree_from_list(mut_pairs):
    parents = []
    children = []
    pairs_of_mutations = []
    for item in mut_pairs:
        a = 'r'+str(item[0])
        b = 'r'+str(item[1])
        pairs_of_mutations.append((a,b))
    t = Tree() # Creates an empty tree
    r0 = t.add_child(name="r0")
    lookup = {"r0": r0}

    for pair in sorted(pairs_of_mutations, key=sort_pairs):
        parentname = pair[0]
        childname = pair[1]
        if childname not in lookup:
            if parentname in lookup:
                newchild = lookup[parentname].add_child(name = childname)
                lookup.update({childname: newchild})

                parents.append(parentname) #make list of unique terminal nodes (no children of children)
                children.append(newchild)
            else:
                print(pair)
                raise RuntimeError('Must not happen.')

    return t

def make_pruned_tree_from_list(mut_pairs):
    parents = []
    children = []
    pairs_of_mutations = []
    for item in mut_pairs:
        a = 'r'+str(item[0])
        b = 'r'+str(item[1])
        pairs_of_mutations.append((a,b))
    t = Tree() # Creates an empty tree
    r0 = t.add_child(name="r0")
    lookup = {"r0": r0}
    prune_list = ['r0']
    for pair in sorted(pairs_of_mutations, key=sort_pairs):
        parentname = pair[0]
        childname = pair[1]
        if childname not in lookup:
            if parentname in lookup:
                newchild = lookup[parentname].add_child(name = childname)
                lookup.update({childname: newchild})
                if parentname not in parents:
                    prune_list.append(lookup[parentname])
                parents.append(parentname) #make list of unique terminal nodes (no children of children)
                children.append(newchild)
            else:
                print(pair)
                raise RuntimeError('Must not happen.')
    prune_count = Counter(children)
    t.prune(prune_list)
    return t


##simulation of bacterial evolution over a rectangle of given size with discrete (horizontal) bands of specified antiobitic concentration
##outputs the cellstate over the entire course of evolution
def megaplate_sim(length ,width, divs,abx,mut_rate, k,s_p,s,pk,density):
    ##set up networkx grid to calculate neighbors
   ##set up networkx grid to calculate neighbors
    G = grid_graph(dim=[width, length])

    ##set up array to track mutation states of cells in grid
    cells = np.full((length,width),-1)

    ##'innoculate' the first row with wild type cells
    cells[0] = 0

    ##actual cells per grid cell, as scaled by the MegaPLate experiment,  Baym et. al 
    cells_per_grid = density

    ##mutation odds over mutation space based on the 
    ##expected maximum value for the number (of cells in grid_cell) drawn from a  poisson distruibtion  
    mut_space = np.linspace(0,pk,pk+1)
    #prob_dist = poisson_max_cdf(np.linspace(0,100,101),mut_rate,cells_per_grid)
    prob_dist = poisson.cdf(np.linspace(0,100,101),mut_rate)
    ##save all cell 
    #cell_history = []

    ##indicate that the half time has been reached yet
    half_time =0 

    count= 0
    
    ##when to stop running the cimulation
    #while all(cells[-1] == -1) and len(cell_history) != 40002:
    #while all(cells[-int(ln/g)+2] == -1) and len(cell_history) != 40002:
    
    #stop running when a cell has made it (two rows in) to the last band
    while all(cells[-int(length/divs)+2] == -1):



        ##find slots where there is a living cells
        cells_where =  np.where(cells >= 0)


        ##create a randomized list of the living cells with which to iterate through
        live_cells_list = []
        for x, y in zip(cells_where[0], cells_where[1]):
            live_cells_list.append([x,y])

        np.random.shuffle(live_cells_list)
        
        
        for j in live_cells_list:
                
                        ##assign some amount (including 0) of additional mutations to the parent cell spot
                        #which now effectively becomes one of the daughter cells

                    
            empty_neighbors = [x for x in G.neighbors(tuple(j)) if cells[tuple((x))] ==-1 ] 
            #if the array containing the empty neighbors is not empty picj one at random
            if empty_neighbors:
                g_draw = 2 * random() -1
                if fitness_func(abx[tuple(j)],cells[tuple(j)],k,s_p,s,pk)> g_draw :
                    pick = choice(empty_neighbors)
                    m = mutation(random(),prob_dist)

                    m_net = []
                    for i in range(m):
                        m_net.append(choice([-1,1]))
                        cells[tuple(j)] = cells[tuple(j)]+sum(m_net)

                    ##other daughter cells occupies chosen empty spot
                    cells[tuple(pick)] = cells[tuple(j)]

                    ##check if half way has been reached, only need to check when a 
                    #division occurs instead of each iteration to save time
                    if all(cells[-3*int(length/divs)+2] == -1) and half_time==0:
                        half_time = count
                    
        count = count +1   
    return count,half_time,cells