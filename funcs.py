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


def fitness_func(d,m,k,s_p,s,pk):
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
def megaplate_sim(length ,width, divs,abx,mut_rate, k,s_p,s,pk,fitness_function):
    ##set up networkx grid to calculate neighbors
    G = grid_graph(dim=[width, length])

    ##set up array to track mutation states of cells in grid
    cells = np.full((length,width),-1)

    ##'innoculate' the first row with wild type cells
    cells[0] = 0


    ##mutation odds over mutation space based on the 
    ##expected maximum value for the number (of cells in grid_cell) drawn from a  poisson distruibtion  
    mut_space = np.linspace(0,pk,pk+1)
    #prob_dist = poisson_max_cdf(np.linspace(0,100,101),mut_rate,cells_per_grid)
    prob_dist = poisson.cdf(np.linspace(0,100,101),mut_rate)


    ##indicate that the half time has been reached yet

    #iniate counter
    count = 0 
    
    #set halfway mark flag to 0
    half_time = 0

    #save all mutation
    mut_count = 0
    mut_states = np.full((length,width),0)
    mut_pairs = []
    
    
    #stop running when a cell has made it to the last band(two gridcells in)
    while all(cells[-int(length/divs)+2] == -1):


        ##find slots where there is a living cells
        cells_where =  np.where(cells >= 0)


        ##create a randomized list of the living cells with which to iterate through
        live_cells_list = []
        for x, y in zip(cells_where[0], cells_where[1]):
            live_cells_list.append([x,y])

        np.random.shuffle(live_cells_list)
        
        
        for j in live_cells_list:
            #cfeate list of tuples of jth cells empty neighbors if any
                    
            empty_neighbors = [x for x in G.neighbors(tuple(j)) if cells[tuple((x))] ==-1 ] 
            #if the array containing the empty neighbors is not empty allow it to mutate
            if empty_neighbors:
                #random number between -1 and 1 
                g_draw = 2 * random() -1
                
                #compare random number to computed fitness to see if fitness is greater than drawn value
                if fitness_function(abx[tuple(j)],cells[tuple(j)],k,s_p,s,pk)> g_draw :
                    #choose an neighbor at random (from list of tuples)
                    pick = choice(empty_neighbors)
                    #number of mutations to be applied to daughter cells
                    m = mutation(random(),prob_dist)*choice([-1,1])
                    
                    if m != 0:
                        #mutations are applied to daughter cell in current grid cell                   
                        cells[tuple(j)] += m
                    
                        ##create list pair of the parent and daughter cells mutation event count
                        mut_pairs.append([mut_states[tuple(j)],mut_count])

                        #update the count of the mutation event that the daughter cell in same grid cell just inherited
                        mut_states[tuple(j)] = mut_count
                        
                        #advance count of mutation event
                        mut_count += 1
                    
                    #update the count of the mutation event that the other daughter cell inherited
                    mut_states[tuple(pick)] = mut_states[tuple(j)]
                                        
                    ##other daughter cells occupies chosen empty spot
                    cells[tuple(pick)] = cells[tuple(j)]

                    ##check if made it to half way
                    if all(cells[-3*int(length/divs)+2]!=-1) and half_time ==0:
                        half_time = count + 1
                    



        count += 1            
    return count,half_time,cells,mut_pairs