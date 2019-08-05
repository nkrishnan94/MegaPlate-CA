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

def sim_user_abx(wt,ln, g, ab, 
                 mut_rate,k, s_p,s,pk,
                 mut_select,track_pairs):
    
    
    
    width = wt
    length = ln
    grid_size = g
    
    ## bacteria per grid cell as mapped to MegaPLate
    bps = int(((60*((120/9)*5) * 10**8)/1)/(ln*wt))
    ##setup nx graph to find neibors quickly
    G = grid_graph(dim=[wt, ln])
    
    ##set up grid to keep track of cell state i.e. number of mutations
    cells = np.full((ln,wt),-1)
    cells[0] = 0 
    if track_pairs == True:
        ##set up grid to keep track of mutation events i.e. which order event resulted in the cell here if any
        muts = np.full((ln,wt),0)
    
    #dstribution of for max of scamples drawn from poisson
    t = np.linspace(0,pk,pk+1)
    if mut_select == 3:
        prob_dist = poisson_max_cdf(t,mut_rate,bps)

    #new 
    if mut_select == 0:

        t = np.linspace(0,pk,pk+1)
        poisson_unbias = poisson.pmf(t,mut_rate)
        poisson_biased=[]
        for i in np.unique(abx_grad):
            poisson_biased.append((poisson.pmf(t,mut_rate)*((fitness(i,t,k,s_p,s,pk)+1)/2))/(sum(poisson.pmf(t,mut_rate)*((fitness(i,t,k,s_p,s,pk)+1)/2))))

    if mut_select == 1:
        t = np.linspace(0,pk,pk+1)
        poisson__max_unbias = poisson_max_cdf(t,mut_rate,bps)
        poisson_biased=[]
        for i in np.unique(abx_grad):
            poisson_biased.append((poisson__max_unbias*((fitness(i,t,k,s_p,s,pk)+1)/2)))   
        
    
    
    
    ##set up grid that maps abx conc.  to space
    ab = ab
    
    #storing cells and mutation
    cell_history = []
    mut_pairs = []
    mut_ID = 0
    half_time = 0
    ##begin evolution
    #while all(cells[-1] == -1) and len(cell_history) != 40002:
    #while all(cells[-int(ln/g)+2] == -1) and len(cell_history) != 40002:
        
    while all(cells[-1] == -1):
        
        ## save current state map to list
        cell_history.append(cells.tolist())
        if half_time ==0:
            if all(cells[-3*int(ln/g)+2] == -1) == False:
                half_time = len(cell_history)
        
        ##find slots where there is a living cells
        cells_where =  np.where(cells != -1)
        
        
        ##create a randomized list of the living cells with which to iterate through
        cells_list = []
        for x, y in zip(cells_where[0], cells_where[1]):
            cells_list.append([x,y])
            
        np.random.shuffle(cells_list)
        
        ##decide if each living in this generation will die, live, or mutate
        for j in cells_list: 
            g_draw = 2 * random() -1
            
            ##death
            j_muts = cells[tuple(j)]
            antibiotic_value = ab[tuple(j)]
            if fitness(antibiotic_value,j_muts,k,s_p,s,pk) < g_draw :
                cells[tuple(j)] == -1
            else:
                neighbors = [x for x in G.neighbors(tuple(j))]

                #find which of the neighboring cells are empty, and divide, with a daughter cell in that space
                empty = np.where(-1 == np.array([cells[tuple(x)] for x in neighbors]) )
                if len(empty[0]) != 0:
                    pick = neighbors[choice(empty)[0]]

                    #mutated daughter cells
                    #m = mutation(random(),prob_dist)
                    #m = np.random.poisson(mut_rate)

                    if mut_select == 2:
                        m = np.random.poisson(mut_rate)
                    if mut_select ==3:
                        m = mutation(random(),prob_dist)
                    if mut_select ==0:
                        m = np.random.choice(t,1,p = poisson_biased[np.where(ab[tuple(j)]== np.unique(ab))[0][0]])

                    if m != 0:
                        if track_pairs == True:
                            mut_ID = mut_ID +1
                            mut_pairs.append([muts[tuple(j)],mut_ID])

                            muts[tuple(j)] = mut_ID
                        cells[tuple(j)] = cells[tuple(j)]+m
                    #divide
                    cells[tuple(pick)] = cells[tuple(j)]
                    


        
    return cell_history, mut_pairs, half_time