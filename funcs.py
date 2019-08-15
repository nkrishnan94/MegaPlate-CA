from networkx import grid_graph, neighbors
import numpy as np
from scipy.stats import poisson 
from random import random, sample, choice

def poisson_max_cdf(x,mu,n):
    y = poisson.cdf(x,mu)**n
    return y

def mutation(r,dist):
    check = np.append(dist,r)
    check.sort()
    m = np.where(check ==r )[0][0]
    return m

def fitness_func(d,m,k,s_p,s,pk):
    p = s - s_p
    y = s + (p)* ((d**k)/((d**k)- (s-p)/s))*((m - pk)/pk)
    return y

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

    count = 0
    
    ##when to stop running the cimulation
    #while all(cells[-1] == -1) and len(cell_history) != 40002:
    #while all(cells[-int(ln/g)+2] == -1) and len(cell_history) != 40002:
    
    #stop running when a cell has made it to the last row
    while all(cells[-1] == -1):
        if all(cells[3* int(length/divs) -2] != -1 and half_time = 0:
               half_time = count 
        ##each itteration the  most recent cell state is added to the cell_history array
        cell_history.append(cells.tolist())




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
                    
            
        count = count+1
    return count,half_time,cells








