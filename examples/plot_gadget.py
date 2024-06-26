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
data=loadtxt("particle_dump.dat")
id = data[:,0]
accx_e=data[:,5]
accy_e=data[:,6]
accz_e=data[:,7]

accx_bh=data[:,8]
accy_bh=data[:,9]
accz_bh=data[:,10]

accx_new=data[:,11]
accy_new=data[:,12]
accz_new=data[:,13]

# Sort accelerations
rank = argsort(id)
id = id[rank]
accx_e = accx_e[rank]
accy_e = accy_e[rank]
accz_e = accz_e[rank]

accx_bh = accx_bh[rank]
accy_bh = accy_bh[rank]
accz_bh = accz_bh[rank]

accx_new = accx_new[rank]
accy_new = accy_new[rank]
accz_new = accz_new[rank]

# Read Gadget accelerations
data=loadtxt("gadget_dump.dat")
id = data[:,0]
accx_g=data[:,5]
accy_g=data[:,6]
accz_g=data[:,7]

# Sort accelerations
rank = argsort(id)
id = id[rank]
accx_g = accx_g[rank]
accy_g = accy_g[rank]
accz_g = accz_g[rank]


# Build error ------------------------------------------------

errx_bh = (accx_bh - accx_e )/abs(accx_e) 
erry_bh = (accy_bh - accy_e )/abs(accy_e) 
errz_bh = (accz_bh - accz_e )/abs(accz_e) 

errx_new = (accx_new - accx_e )/abs(accx_e) 
erry_new = (accy_new - accy_e )/abs(accy_e) 
errz_new = (accz_new - accz_e )/abs(accz_e) 

errx_g = (accx_g - accx_e )/abs(accx_e) 
erry_g = (accy_g - accy_e )/abs(accy_e) 
errz_g = (accz_g - accz_e )/abs(accz_e) 

# Statistics
meanx_bh = mean(errx_bh[abs(errx_bh) < 0.1])
stdx_bh = std(errx_bh[abs(errx_bh) < 0.1])
meany_bh = mean(erry_bh[abs(erry_bh) < 0.1])
stdy_bh = std(erry_bh[abs(erry_bh) < 0.1])
meanz_bh = mean(errz_bh[abs(errz_bh) < 0.1])
stdz_bh = std(errz_bh[abs(errz_bh) < 0.1])

meanx_new = mean(errx_new[abs(errx_new) < 0.1])
stdx_new = std(errx_new[abs(errx_new) < 0.1])
meany_new = mean(erry_new[abs(erry_new) < 0.1])
stdy_new = std(erry_new[abs(erry_new) < 0.1])
meanz_new = mean(errz_new[abs(errz_new) < 0.1])
stdz_new = std(errz_new[abs(errz_new) < 0.1])

# meanx_g = mean(errx_g)
# stdx_g = std(errx_g)
# meany_g = mean(erry_g)
# stdy_g = std(erry_g)
# meanz_g = mean(errz_g)
# stdz_g = std(errz_g)

# Plot -------------------------------------------------------
figure(frameon=True)

subplot(311, title="Acceleration along X")
plot(id, errx_g , 'gs')
plot(id, errx_bh , 'rx')
plot(id, errx_new , 'b.')
text(id[-1], 0.18, "B-H: $%5.3f\\pm%5.3f$\n QuickShed: $%5.3f\\pm%5.3f$"%(meanx_bh, stdx_bh, meanx_new, stdx_new), backgroundcolor="w", va="top", ha="right" )


ylim(-0.2, 0.2)
xlim(0,id[-1])
grid()

subplot(312, title="Acceleration along Y")
plot(id, erry_g , 'gs')
plot(id, erry_bh , 'rx')
plot(id, erry_new , 'b.')
text(id[-1], 0.18, "B-H: $%5.3f\\pm%5.3f$\n QuickShed: $%5.3f\\pm%5.3f$"%(meany_bh, stdy_bh, meany_new, stdy_new), backgroundcolor="w", va="top", ha="right" )


ylim(-0.2, 0.2)
xlim(0,id[-1])

grid()

subplot(313, title="Acceleration along Z")
plot(id, errz_g , 'gs', label="Gadget")
plot(id, errz_new , 'b.', label="QuickShed")
plot(id, errz_bh , 'rx', label="Legacy")
#text(id[-1], 0.18, "B-H: $%5.3f\\pm%5.3f$\n QuickShed: $%5.3f\\pm%5.3f$"%(meanz_bh, stdz_bh, meanz_new, stdz_new), backgroundcolor="w", va="top", ha="right" )
legend(loc="upper right")

ylim(-0.2, 0.2)
xlim(0,id[-1])
grid()

savefig("accelerations.png")




# Plot -------------------------------------------------------
bins = linspace(-3, 3, 10000)


figure(frameon=True)
subplot(311, title="Acceleration along X")#, yscale='log')
hist(errx_g, bins=bins, normed=1, histtype='step', rwidth=0.01, color='g', label="Gadget")
hist(errx_bh, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r', label="Legacy")
hist(errx_new, bins=bins, normed=1, histtype='step', rwidth=0.01, color='b', label="QuickShed")
legend(loc="upper right")
xlim(-0.03, 0.03)

subplot(312, title="Acceleration along Y")
hist(erry_g, bins=bins, normed=1, histtype='step', rwidth=0.01, color='g')
hist(erry_bh, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r')
hist(erry_new, bins=bins, normed=1, histtype='step', rwidth=0.01, color='b')
xlim(-0.03, 0.03)

subplot(313, title="Acceleration along Z")
hist(errz_g, bins=bins, normed=1, histtype='step', rwidth=0.01, color='g')
hist(errz_bh, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r')
hist(errz_new, bins=bins, normed=1, histtype='step', rwidth=0.01, color='b')
xlim(-0.03, 0.03)

savefig("histogram.png")


# Print error
errx_g = errx_g[ (errx_g > -0.1) & (errx_g < 0.1) ]
erry_g = erry_g[ (erry_g > -0.1) & (erry_g < 0.1) ]
errz_g = errz_g[ (errz_g > -0.1) & (errz_g < 0.1) ]

errx_bh = errx_bh[ (errx_bh > -0.1) & (errx_bh < 0.1) ]
erry_bh = erry_bh[ (erry_bh > -0.1) & (erry_bh < 0.1) ]
errz_bh = errz_bh[ (errz_bh > -0.1) & (errz_bh < 0.1) ]

errx_new = errx_new[ (errx_new > -0.1) & (errx_new < 0.1) ]
erry_new = erry_new[ (erry_new > -0.1) & (erry_new < 0.1) ]
errz_new = errz_new[ (errz_new > -0.1) & (errz_new < 0.1) ]


print "Error for Gadget: ( x= %4.3e"%std(errx_g), "y= %4.3e"%std(erry_g), "z= %4.3e"%std(errz_g), ") Combined= %4.3e"%(( std(errx_g) + std(erry_g) + std(errz_g) )/3.)
print "Error for Legacy: ( x= %4.3e"%std(errx_bh), "y= %4.3e"%std(erry_bh), "z= %4.3e"%std(errz_bh), ") Combined= %4.3e"%(( std(errx_bh) + std(erry_bh) + std(errz_bh) )/3.)
print "Error for QuShed: ( x= %4.3e"%std(errx_new), "y= %4.3e"%std(erry_new), "z= %4.3e"%std(errz_new), ") Combined= %4.3e"%(( std(errx_new) + std(erry_new) + std(errz_new) )/3.)
