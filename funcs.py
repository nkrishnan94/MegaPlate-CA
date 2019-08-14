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

def megaplate_sim(length,width, divs,abx,mut_rate, k,s_p,s,pk):
    ##set up networkx grid to calculate neighbors
   ##set up networkx grid to calculate neighbors
    G = grid_graph(dim=[width, length])

    ##set up array gto track mutation states of cells in grid
    cells = np.full((length,width),-1)

    ##'innoculate' the first row with wild type cells
    cells[0] = 0

    ##actual cells per grid cell, as scaled by the MegaPLate experiment,  Baym et. al 
    cells_per_grid = int(((60*((120/9)*5) * 10**8)/1)/(length*width))

    ##mutation odds over mutation space based on the 
    ##expected maximum value for the number (of cells in grid_cell) drawn from a  poisson distruibtion  
    mut_space = np.linspace(0,pk,pk+1)
    prob_dist = poisson_max_cdf(np.linspace(0,100,101),mut_rate,cells_per_grid)

    ##save all cell 
    #cell_history = []

    ##indicate that the half time has been reached yet
    half_time =0 

    cell_history= []
    count = 0 
    ##when to stop running the cimulation
    #while all(cells[-1] == -1) and len(cell_history) != 40002:
    #while all(cells[-int(ln/g)+2] == -1) and len(cell_history) != 40002:

    while all(cells[-1] == -1):
        ##each itteration the  most recent cell state is added to the cell_history array
        cell_history.append(cells.tolist())
        #count = count + 1
        ##check if half way mark is reached 
        if half_time ==0:
            if all(cells[-3*int(length/divs)+2] == -1) == False:
                half_time = count



        ##find slots where there is a living cells
        cells_where =  np.where(cells != -1)


        ##create a randomized list of the living cells with which to iterate through
        cells_list = []
        for x, y in zip(cells_where[0], cells_where[1]):
            cells_list.append([x,y])

        np.random.shuffle(cells_list)


        for j in cells_list: 
            g_draw = 2 * random() -1

            ##death
            #fitness = fitness_values[fitness_values['mutations']==cells[tuple(j)]][fitness_values['drug']==ab[tuple(j)]]['fitness'].iloc[0]
            if fitness_func(abx[tuple(j)],cells[tuple(j)],k,s_p,s,pk)< g_draw :
                cells[tuple(j)] == -1
            else:
                empty_neighbors = [x for x in G.neighbors(tuple(j)) if cells[tuple((x))] ==-1 ] 
                if empty_neighbors: 
                    pick = choice(empty_neighbors)


                    m = mutation(random(),prob_dist)

                    if m!=0:
                        cells[tuple(j)] = cells[tuple(j)]+m


                    cells[tuple(pick)] = cells[tuple(j)]


        count = count + 1
    return cell_history, half_time, cells








