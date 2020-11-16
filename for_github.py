import numpy as np
from numpy import genfromtxt
from matplotlib.pylab import *
import matplotlib.pyplot as plt
import os

#presa dati
file=input('file name? ')
data = genfromtxt(file, delimiter='	', skip_header=1, usecols=(1,2,3,4,5,6))
R, z, phi, vR, vz, vphi = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5]

n_data=R.size


#funzioni

def dist2( _R,_vR, _i,_z=np.zeros(n_data), _phi=np.zeros(n_data),  _vz=np.zeros(n_data), _vphi=np.zeros(n_data)): 
    w=(sqrt((_R-_R[_i])**2+(_z-_z[_i])**2+(_phi-_phi[_i])**2+(_vR-_vR[_i])**2+(_vz-_vz[_i])**2+(_vphi-_vphi[_i])**2))
    return w

def get_param(_x,_y):
    N=len(_x)
    #M=number of parameters
    sigma=np.zeros(N)

    #for i in range(N):
    #sigma[i]=1 	in generale (y_i-y(x_i))**/(N-M)

    sigma=np.ones(N)

    S=sum(1/sigma**2)
    S_x=sum(_x/sigma**2)
    S_y=sum(_y/sigma**2)
    S_xx=sum(_x**2/sigma**2)
    S_xy=sum(_x*_y/sigma**2)
    delta=S*S_xx - S_x**2

    #parameters
    q=(S_xx*S_y - S_x*S_xy)/delta
    eq=np.sqrt(S_xx/delta)

    #parameters' errors
    m=(S*S_xy - S_x*S_y)/delta
    em=np.sqrt(S/delta)

    return np.array([m,q]), np.array([em,eq])

def min_chi_square(_x, _y):
    #creating the window
    L=len(_x)
    l=500
    step = L - l

    #emptyarrays
    param=np.zeros((step,2))
    errparam=np.zeros((step,2))
    chi_square=np.zeros(step)

    for i in range(step):
        X=_x[i:l+i]
        Y=_y[i:l+i]
        param[i],errparam[i]=get_param(X,Y)
        chi_square[i]=sum((Y-param[i][0]*X -param[i][1])**2)

    index_abs_minimum = np.argmin(chi_square)

    index_first_minimum=0
    while chi_square[index_first_minimum]>chi_square[index_first_minimum+1]:
        index_first_minimum+=1

    #adjusting shape
    parameters=np.ravel(param[index_abs_minimum])
    errparameters=np.ravel(errparam[index_abs_minimum])
    bestparameters=np.ravel(param[index_first_minimum])
    errbestparameters=np.ravel(errparam[index_first_minimum])

    #output
    print("Best coeff ang :", bestparameters[0], '+/-', errbestparameters[0])
    print("Best normalization :", bestparameters[1], '+/-', errbestparameters[1])

    #writing chi square data on txt file
    chidata = np.array([_x[:step], chi_square, param[:,0]])
    chidata = chidata.T

    chipath="chisquare"

    if not os.path.exists(chipath):
        os.mkdir(chipath)

    chifile="chidata"+file+".txt"
    datafile_path = "chisquare/"+chifile
    with open(datafile_path, 'w') as datafile_id:
        np.savetxt(datafile_id, chidata, delimiter='    ', header='r    chi    dim', comments='')

    return parameters, bestparameters, errbestparameters


#computer distances between all points
distances=np.zeros((n_data,n_data)) 

for i in range(n_data):
    distances[i]= dist2(R,vR,i)

#creating the histogram
distances = distances.reshape(n_data*n_data)
grid=np.logspace(-6, np.log10(np.amax(distances)),int(n_data/2))
results = np.histogram(distances, bins=grid)

#compute C(r) with cumulative histogram
rr = (results[1][1:] + results[1][:-1])/2.
cumulative = np.zeros(len(rr))
cumulative[0] = results[0][0]

for i in range(1,len(rr),1):
    cumulative[i] = cumulative[i-1]+results[0][i]

cumulative = cumulative/n_data**2

#fit
s=cumulative[np.where(cumulative>5*10**(-5))]
r=rr[np.where(cumulative>5*10**(-5))]

param, bestparam, errbestparam=min_chi_square(log10(r),log10(s))

m=np.zeros(rr.size)
m[0]= bestparam[0]
q=np.zeros(rr.size)
q[0]= bestparam[1]

#writing correlation integral data on file
corr_int_data = np.array([rr, cumulative, m, q])
corr_int_data = corr_int_data.T

corr_int_path="correlation_integral"
if not os.path.exists(corr_int_path):
        os.mkdir(corr_int_path)

corr_int_file="corr_int_data"+file+".txt"
corr_datafile_path = "correlation_integral/"+corr_int_file

with open(corr_datafile_path, 'w') as corr_datafile_id:
    np.savetxt(corr_datafile_id, corr_int_data, delimiter='    ', header='r    corr_int    m    q', comments='')
