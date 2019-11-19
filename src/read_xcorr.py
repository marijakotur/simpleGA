'''
Created on 19 sep. 2016

@author: Marija
'''
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

dir_name = 'C:/work/2016_xcorr'
file_name = 'scandata_2016-09-09_14h37.txt'
whole_path_name = dir_name + '/' + file_name

time_data = []
mm_data = []
counts_data = []
datafile = open(whole_path_name)
for line in datafile:
    a = line.split()
    mm_data.append(float(a[0]))
    time_data.append(float(a[1])*1e12)
    counts_data.append(float(a[2]))

plt.subplot()
plt.plot(time_data,counts_data)

plt.xlabel('Time [ps]')
plt.ylabel('PD signal [arb]')



#fitting polynomial
f = np.poly1d(np.polyfit(time_data, counts_data, 16))
fitted = f(time_data)
plt.plot(time_data,fitted)

#fit arbitrary function


def gfun(x, *args):
    #return a * np.exp(-(x-b)**2/(2*c**2)) + d
    c = []
    for a in args[0]:
        c.append(a)
    print c
    #return c[0] * np.exp(-(x-c[1])**2/(2*c[2]**2)) + c[3]

def gfun4(x, a, w, b1, b2,b3, b4, c):
    return a * np.exp(-(x-b1)**2/(2*w**2)) + np.exp(-(x-b2)**2/(2*w**2)) + np.exp(-(x-b3)**2/(2*w**2)) + np.exp(-(x-b4)**2/(2*w**2)) + c

a = np.array([1, 7, 2.5, 0.8])
initial = gfun(time_data,a)
#plt.plot(time_data,initial)

popt, pcov = curve_fit(gfun4, time_data, counts_data, [0.25, 1.5, np.mean(time_data)-1.5,np.mean(time_data)-0.5,np.mean(time_data)+0.5,np.mean(time_data)+1.5,1])
print popt
fitted = gfun4(time_data, popt[0], popt[1],popt[2],popt[3],popt[4],popt[5],popt[6])

plt.plot(time_data,fitted)

plt.show()
