from funcs import sim_user_abx
import numpy as np 
import shelve
import tqdm


ln = [80,140]
wd = [int(.9*i) for i in ln]

gr = 5

<<<<<<< Updated upstream
=======
length = [40, 60]
log_file.write("length: [" )
for item in length:
	log_file.write("%f, " % item)
log_file.write("]\n" )
>>>>>>> Stashed changes




mut_rate = [.001]

k = [4,8]
s_p= -1
s = 1
pk = [200]

<<<<<<< Updated upstream
=======
mut_rate = [.001,.01]
log_file.write("mutations rates: [" )
for item in mut_rate:
	log_file.write("%f, " % item)
log_file.write("]\n" )
>>>>>>> Stashed changes

reps = 1
MegaPlate_42919 = []
    
for a in range(reps): 
    for b in range(len(mut_rate)):
        for c in range(len(ln)):
            for d in range(len(k)):
                for e in tqdm.tqdm(range(len(pk))):
                    
                    l = [val for sublist in [[j for i in range(int(ln[c]/gr))] for j in range (0,gr)] for val in sublist]
                    ab =[[3*10**(i-1)]*wd[c] for i  in l if i != 0]
                    [ab.insert(i, [0]*wd[c]) for i in range(0,int(ln[c]/gr))]
                    abx_grad = np.array(ab)

                



                    cellhistory,mut_pairs,half_time =sim_user_abx(wt = int(wd[c]),ln = int(ln[c]), g = gr,ab = abx_grad, 
                                                        mut_rate= mut_rate[b],k = k[d], s_p = s_p,
                                                        s = s, pk=pk[e],mut_select = 0,track_pairs=False)

                    time = (len(cellhistory))
                    muts = np.array(cellhistory).max()
                    MegaPlate_72919.append({"reps":a,"mut_rate":mut_rate[b],"length":ln[c],
                                                        "k":k[d],'pk':pk[e],"mut_select": 0,"time":time,
                                                        "muts":muts,"half_time":half_time})
                    shelf = shelve.open("MegaPlate_72919.shlf")

                    # serializing
                    shelf["my_dict"] = MegaPlate_exhausitive_42419

                    shelf.close()

<<<<<<< Updated upstream

=======
reps = 1
log_file.write("number of repitions for each set of parameters: %f\n" % reps)

total_runs = reps*len(pk)*len(k)*len(mut_rate)*len(length)*len(width)
MegaPlate_7619 = []




i = 0     
for a in range(reps): 
    for b in range(len(length)):
    	for c in range(len(width)):
        	for d in range(len(mut_rate)):
        		for e in range(len(k)):
        			for f in range(len(pk)):

        				l = [val for sublist in [[j for i in range(int(length[b]/divs))] for j in range (0,divs)] for val in sublist]
        				ab =[[3*10**(i-1)]*width[c] for i  in l if i != 0]
        				[ab.insert(i, [0]*width[c]) for i in range(0,int(length[b]/divs))]
        				abx_grad = np.array(ab)


        				t,t_half,cells =megaplate_sim(width = int(width[c]),length = int(length[b]), divs=divs,abx = abx_grad,
			                                                          mut_rate= mut_rate[b],k = k[e], s_p = s_p,
			                                                          s = s, pk=pk[f])
        				time = t
        				muts = np.array(cells).max()
        				log_file.write('%d runs out of %d completed\n' % (i,total_runs))
        				i = i+1

                        MegaPlate_7619.append({"reps":a,"mut_rate":mut_rate[d],"length":length[b],"k":k[e],'pk':pk[f],"time":t,"muts":muts,"half_time":t_half})


# file to be useds
shlf_file = "MegaPlate_7619.shlf"
shelf = shelve.open("%s" % shlf_file)

# serializing
shelf["my_dict"] = MegaPlate_7619
shelf.close()
log_file.write("results saved: %s" % shlf_file)
log_file.close()
>>>>>>> Stashed changes


