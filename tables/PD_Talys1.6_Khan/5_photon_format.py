from numpy import *
import os

# -----------------------------------------------
# Renormalize the cross sections of the photon emission channels which remain after thinning
# that their sum before and after the thinning is equal.
# Constrain cross section values to inelastic cross section values if cross section is larger
# than inelastic cross section
# (i.e. Number of photo disintegrations from mother to excited daughter nucleus with emission 
# of photon of specific energy can not exceed number of photo disintegrations from mother to
# daughter nucleus).
# -----------------------------------------------

data = genfromtxt('xs_excl_photon.txt')
daughter = vstack({tuple(r) for r in data[:,:4]})
daughter = daughter.astype(int)
daughter = sort(daughter.view('i8,i8,i8,i8'), order=['f0','f1','f2','f3'],axis=0).view(int)
data_thin = genfromtxt('xs_photon_thin.txt')
data_all_daughters = genfromtxt('xs_sum_daughter.txt')

for (Z,N,Zd,Nd) in daughter:
    print Z, N, Zd, Nd
    idx = (data[:,0] == Z) * (data[:,1] == N) * (data[:,2] == Zd) * (data[:,3] == Nd)
    idy = (data_thin[:,0] == Z) * (data_thin[:,1] == N) * (data_thin[:,2] == Zd) * (data_thin[:,3] == Nd)
    xs = data[:,5:][idx]
    xs_thin = data_thin[:,5:][idy]
    xs_sum = sum(xs, axis=0)
    xs_sum_thin = sum(xs_thin, axis=0)

    # renormalize cross sections
    renormalization = xs_sum/xs_sum_thin
    renormalization[isnan(renormalization)] = 0
    data_thin[:,5:][idy] *= renormalization 

    # constrain cross sections
    idz = (data_all_daughters[:,0] == Z) * (data_all_daughters[:,1] == N) * (data_all_daughters[:,2] == Zd) * (data_all_daughters[:,3] == Nd)
    xs = data_all_daughters[:,4:][idz]
    xs = xs[0,:]

    temp = data_thin[idy]
    for i in range(0, len(temp[:,0])):
        for j in range(0, len(temp[0,5:])):
            if (temp[i,5+j] > xs[j]):
                temp[i,5+j] = xs[j]
    data_thin[idy] = temp

savetxt('xs_photon_norm.txt', data_thin, fmt='%i\t%i\t%i\t%i\t%.4f' + '\t%g'*301)
