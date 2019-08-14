from funcs import sim_user_abx
import matplotlib.pyplot as plt
import numpy as np

ln = 140
wd = int(.9 *ln)
gr = 5

l = [val for sublist in [[j for i in range(int(ln/gr))] for j in range (0,gr)] for val in sublist]
ab =[[3*10**(i-1)]*wd for i  in l if i != 0]
[ab.insert(i, [0]*wd) for i in range(0,int(ln/gr))]
abx_grad = np.array(ab)

ch1,mp,ht = sim_user_abx(wt = wd,ln =ln,g= gr, 
	ab = abx_grad ,mut_rate=.01,
	k = 3, s_p = -1,s = 1, 
	pk= 1000,mut_select = 3
     ,track_pairs = False)

#print(len(ch1))
#plt.imshow(ch1[-1])
#plt.show()

