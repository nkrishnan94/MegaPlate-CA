from funcs import megaplate_sim
import numpy as np
import shelve
import datetime


log_file = open("log_%s.txt" % datetime.datetime.now(),"w+")
log_file.write("%s\r\n" % datetime.datetime.now())


length = [40, 50]
log_file.write("length: [" )
for item in length:
    log_file.write("%f, " % item)
log_file.write("]\n" )



width = [int(.9*i) for i in length]
log_file.write("width: [" )
for item in width:
    log_file.write("%f, " % item)
log_file.write("]\n" )

divs = 5
log_file.write("Antibiotic divisions: %f\n" % divs )


mut_rate = [.001,.01]
log_file.write("mutations rates: [" )
for item in mut_rate:
    log_file.write("%f, " % item)
log_file.write("]\n" )

k = [0.4,0.8]

log_file.write("k, or 'hill-like' constant: [" )
for item in k:
    log_file.write("%f, " % item)
log_file.write("]\n" )


s_p= -1
log_file.write("s_p, or minimum fitness: %f\n" % s_p )
s = 1
log_file.write("s_p, or maximum fitness: %f\n" % s )
pk = [200,300]

log_file.write("pk, or maximum possible mutations: [" % pk)
for item in pk:
    log_file.write("%f, " % item)
log_file.write("]\n" )

reps = 1
log_file.write("number of repitions fo each set of parameters: %f\n" % reps)

density = 

tot_runs = reps*len(length)*len(width)*len(mut_rate)*len(k)*len(pk)
##empty list for each dictiorary result to be saved as db file
MegaPlate_7619 = []

#counter for log files
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
                                                                      s = s, pk=pk[f],density)
                        time = t
                        muts = np.array(cells).max()
                        log_file.write('%d out of %d runs completed \n' % (i,tot_runs))
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


