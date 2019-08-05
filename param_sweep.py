from funcs import sim_user_abx
import numpy as np 
import shelve
import tqdm


ln = [80,140]
wd = [int(.9*i) for i in ln]

gr = 5





mut_rate = [.001]

k = [4,8]
s_p= -1
s = 1
pk = [200]


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




