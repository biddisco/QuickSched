#!/usr/bin/env python2                                                                                                                     
# -*- coding: utf-8 -*-                                                                                                                    
import matplotlib
matplotlib.use("Agg")
from pylab import *
import numpy
# tex stuff                                                                                                                                
#rc('font',**{'family':'serif','serif':['Palatino']})                                                                                      
params = {'axes.labelsize': 16,
          'axes.titlesize': 20,
          'text.fontsize': 14,
          'legend.fontsize': 14,
          'xtick.labelsize': 14,
          'ytick.labelsize': 14,
          'text.usetex': True,

          'figure.figsize' : (12,9),
          'figure.subplot.left'    : 0.07,  # the left side of the subplots of the figure                                                  
          'figure.subplot.right'   : 0.98  ,  # the right side of the subplots of the figure                                               
          'figure.subplot.bottom'  : 0.09  ,  # the bottom of the subplots of the figure                                                   
          'figure.subplot.top'     : 0.9  ,  # the top of the subplots of the figure                                                       
          'figure.subplot.wspace'  : 0.26  ,  # the amount of width reserved for blank space between subplots                              
          'figure.subplot.hspace'  : 0.22  ,  # the amount of height reserved for white space between subplots                             

          'lines.markersize' : 6,
          'lines.linewidth' : 2,

#          'axes.formatter.limits' : (-2, 1),                                                                                              

          'text.latex.unicode': True
        }
rcParams.update(params)
rc('font', family='serif')
import sys
import os

print "Plotting..."

# Read Quickshed accelerations
data=loadtxt("interaction_dump.dat")
id = data[:,0]
accx_u=data[:,5]
accy_u=data[:,6]
accz_u=data[:,7]

accx_s=data[:,8]
accy_s=data[:,9]
accz_s=data[:,10]




# Build error ------------------------------------------------

errx_s = (accx_s - accx_u )/abs(accx_u) 
erry_s = (accy_s - accy_u )/abs(accy_u) 
errz_s = (accz_s - accz_u )/abs(accz_u) 

e_errx_s = errx_s#[abs(errx_s) > 0.001]
e_erry_s = erry_s#[abs(erry_s) > 0.001]
e_errz_s = errz_s#[abs(errz_s) > 0.001]

# Statistics
meanx_s = mean(errx_s[abs(errx_s) < 0.1])
stdx_s = std(errx_s[abs(errx_s) < 0.1])
meany_s = mean(erry_s[abs(erry_s) < 0.1])
stdy_s = std(erry_s[abs(erry_s) < 0.1])
meanz_s = mean(errz_s[abs(errz_s) < 0.1])
stdz_s = std(errz_s[abs(errz_s) < 0.1])



# Plot -------------------------------------------------------
figure(frameon=True)

subplot(311, title="Acceleration along X")
#plot(id[abs(errx_s) > 0.001], e_errx_s , 'rx')
plot(id, e_errx_s , 'rx')
#text(id[-1], 0.18, "B-H: $%5.3f\\pm%5.3f$\n QuickShed: $%5.3f\\pm%5.3f$"%(meanx_bh, stdx_bh, meanx_new, stdx_new), backgroundcolor="w", va="top", ha="right" )
ylim(-0.2, 0.2)
xlim(0, size(id)-1)
grid()

subplot(312, title="Acceleration along Y")
#plot(id[abs(erry_s) > 0.001], e_erry_s , 'rx')
plot(id, e_erry_s , 'rx')
#text(id[-1], 0.18, "B-H: $%5.3f\\pm%5.3f$\n QuickShed: $%5.3f\\pm%5.3f$"%(meany_bh, stdy_bh, meany_new, stdy_new), backgroundcolor="w", va="top", ha="right" )
ylim(-0.2, 0.2)
xlim(0, size(id)-1)

grid()

subplot(313, title="Acceleration along Z")
#plot(id[abs(errz_s) > 0.001], e_errz_s , 'rx', label="Sorted")
plot(id, e_errz_s , 'rx', label="Sorted")
#text(id[-1], 0.18, "B-H: $%5.3f\\pm%5.3f$\n QuickShed: $%5.3f\\pm%5.3f$"%(meanz_bh, stdz_bh, meanz_new, stdz_new), backgroundcolor="w", va="top", ha="right" )
legend(loc="upper right")

ylim(-0.2, 0.2)
xlim(0, size(id)-1)
grid()

savefig("accelerations.png")




# Plot -------------------------------------------------------
bins = linspace(-3, 3, 10000)

e_errx_s = errx_s[(abs(errx_s) > 0.001) & (abs(errx_s) < 1.)]
e_erry_s = erry_s[(abs(erry_s) > 0.001) & (abs(erry_s) < 1.)]
e_errz_s = errz_s[(abs(errz_s) > 0.001) & (abs(errz_s) < 1.)]



figure(frameon=True)
subplot(311, title="Acceleration along X")#, yscale='log')
hist(e_errx_s, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r', label="Sorted")
legend(loc="upper right")
xlim(-0.12, 0.12)

subplot(312, title="Acceleration along Y")
hist(e_erry_s, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r')
xlim(-0.12, 0.12)

subplot(313, title="Acceleration along Z")
hist(e_errz_s, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r')
xlim(-0.12, 0.12)

savefig("histogram.png")





print "E_error for sorted: ( x= %4.3e"%std(e_errx_s), "y= %4.3e"%std(e_erry_s), "z= %4.3e"%std(e_errz_s), ") Combined YZ = %4.3e"%(( std(e_erry_s) + std(e_errz_s) )/2.)

#print std(e_errx_s), std(e_erry_s), std(e_errz_s)


#print "Min   for sorted: ( x= %4.3e"%min(errx_s), "y= %4.3e"%min(erry_s), "z= %4.3e"%min(errz_s), ") Combined YZ = %4.3e"%(min( min(erry_s), min(errz_s) ))
#print "Max   for sorted: ( x= %4.3e"%max(errx_s), "y= %4.3e"%max(erry_s), "z= %4.3e"%max(errz_s), ") Combined YZ = %4.3e"%(max( max(erry_s), max(errz_s) ))



