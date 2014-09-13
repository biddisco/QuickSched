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
plot(id, errx_s , 'rx')
#text(id[-1], 0.18, "B-H: $%5.3f\\pm%5.3f$\n QuickShed: $%5.3f\\pm%5.3f$"%(meanx_bh, stdx_bh, meanx_new, stdx_new), backgroundcolor="w", va="top", ha="right" )
ylim(-0.2, 0.2)
xlim(0, size(id)-1)
grid()

subplot(312, title="Acceleration along Y")
plot(id, erry_s , 'rx')
#text(id[-1], 0.18, "B-H: $%5.3f\\pm%5.3f$\n QuickShed: $%5.3f\\pm%5.3f$"%(meany_bh, stdy_bh, meany_new, stdy_new), backgroundcolor="w", va="top", ha="right" )
ylim(-0.2, 0.2)
xlim(0, size(id)-1)

grid()

subplot(313, title="Acceleration along Z")
plot(id, errz_s , 'rx', label="QuickShed")
#text(id[-1], 0.18, "B-H: $%5.3f\\pm%5.3f$\n QuickShed: $%5.3f\\pm%5.3f$"%(meanz_bh, stdz_bh, meanz_new, stdz_new), backgroundcolor="w", va="top", ha="right" )
legend(loc="upper right")

ylim(-0.2, 0.2)
xlim(0, size(id)-1)
grid()

savefig("accelerations.png")




# Plot -------------------------------------------------------
# bins = linspace(-3, 3, 10000)


# figure(frameon=True)
# subplot(311, title="Acceleration along X")#, yscale='log')
# hist(errx_s, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r', label="Legacy")
# legend(loc="upper right")
# xlim(-0.03, 0.03)

# subplot(312, title="Acceleration along Y")
# hist(erry_s, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r')
# xlim(-0.03, 0.03)

# subplot(313, title="Acceleration along Z")
# hist(errz_s, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r')
# xlim(-0.03, 0.03)

# savefig("histogram.png")



