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

          'lines.markersize' : 2,
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
data=loadtxt( "interaction_dump_all.dat" )
id = data[:,0]

mass = data[:,1]
pos_x = data[:,2]
pos_y = data[:,3]
pos_z = data[:,4]

accx_u=data[:,5]
accy_u=data[:,6]
accz_u=data[:,7]

accx_s=data[:,8]
accy_s=data[:,9]
accz_s=data[:,10]

accx_bh=data[:,11]
accy_bh=data[:,12]
accz_bh=data[:,13]
    

# Build error ------------------------------------------------
inv_acc_tot = 1.0 / sqrt(accx_u**2 + accy_u**2 + accz_u**2)     
errx_s = (accx_s - accx_u ) * inv_acc_tot
erry_s = (accy_s - accy_u ) * inv_acc_tot
errz_s = (accz_s - accz_u ) * inv_acc_tot

errx_bh = (accx_bh - accx_u ) * inv_acc_tot
erry_bh = (accy_bh - accy_u ) * inv_acc_tot
errz_bh = (accz_bh - accz_u ) * inv_acc_tot

e_errx_s = errx_s#[abs(errx_s) > 0.001]
e_erry_s = erry_s#[abs(erry_s) > 0.001]
e_errz_s = errz_s#[abs(errz_s) > 0.001]

e_errx_bh = errx_bh#[abs(errx_bh) > 0.001]
e_erry_bh = erry_bh#[abs(erry_bh) > 0.001]
e_errz_bh = errz_bh#[abs(errz_bh) > 0.001]


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
#plot(id[abs(errx_s) > 0.001], e_errx_s , 'ro')
plot(pos_x, e_errx_bh , 'bx', label='B-H')
plot(pos_x, e_errx_s , 'ro')
#text( 0., 0.1, "axis=( %d %d %d )"%(axis[orientation*3 + 0], axis[orientation*3 + 1], axis[orientation*3 + 2]) , ha='center', backgroundcolor='w', fontsize=14)
ylim(-0.1, 0.1)
grid()

subplot(312, title="Acceleration along Y")
#plot(id[abs(erry_s) > 0.001], e_erry_s , 'ro')
plot(pos_x, e_erry_bh , 'bx', label='B-H')
plot(pos_y, e_erry_s , 'ro')
#text( 0., 0.1, "axis=( %d %d %d )"%(axis[orientation*3 + 0], axis[orientation*3 + 1], axis[orientation*3 + 2]) , ha='center', backgroundcolor='w', fontsize=14)
ylim(-0.1, 0.1)  
grid()

subplot(313, title="Acceleration along Z")
#plot(id[abs(errz_s) > 0.001], e_errz_s , 'ro', label="Sorted")
plot(pos_x, e_errz_bh , 'bx', label='B-H')
plot(pos_z, e_errz_s , 'ro', label="Sorted")
#text( 0., 0.1, "axis=( %d %d %d )"%(axis[orientation*3 + 0], axis[orientation*3 + 1], axis[orientation*3 + 2]) , ha='center', backgroundcolor='w', fontsize=14)

legend(loc="upper right")

ylim(-0.1, 0.1)
#xlim(0, size(id)-1)
grid()

savefig("accelerations_relative_all.png")
close()


# Plot -------------------------------------------------------
figure(frameon=True)

subplot(311, title="Acceleration along X")
#plot(id[abs(errx_s) > 0.001], e_errx_s , 'ro')
plot(pos_x, accx_bh, 'bx')
plot(pos_x, accx_s , 'ro')
#text( 0., 0.1, "axis=( %d %d %d )"%(axis[orientation*3 + 0], axis[orientation*3 + 1], axis[orientation*3 + 2]) , ha='center', backgroundcolor='w', fontsize=14)
ylim(-250, 250)
grid()

subplot(312, title="Acceleration along Y")
#plot(id[abs(erry_s) > 0.001], e_erry_s , 'ro')
#plot(pos_y, accy_u , 'gs')
plot(pos_y, accy_s , 'ro')
plot(pos_y, accz_bh , 'bx', label="B-H")
#text( 0., 0.1, "axis=( %d %d %d )"%(axis[orientation*3 + 0], axis[orientation*3 + 1], axis[orientation*3 + 2]) , ha='center', backgroundcolor='w', fontsize=14)
#ylim(-70, 70)
ylim(-250, 250)
grid()

subplot(313, title="Acceleration along Z")
#plot(id[abs(errz_s) > 0.001], e_errz_s , 'ro', label="Sorted")
plot(pos_z, accz_bh , 'bx', label="B-H")
plot(pos_z, accz_s , 'ro', label="Sorted")
#plot(pos_z, accz_u , 'gs', label="Unsorted")

#text( 0., 0.1, "axis=( %d %d %d )"%(axis[orientation*3 + 0], axis[orientation*3 + 1], axis[orientation*3 + 2]) , ha='center', backgroundcolor='w', fontsize=14)

legend(loc="upper right")

#ylim(-70, 70)
ylim(-250, 250)
grid()

savefig("accelerations_absolute_all.png")
close()



# Plot -------------------------------------------------------
bins = linspace(-3, 3, 10000)

e_errx_s1 = errx_s[ pos_x < 0.5 ]
e_erry_s1 = erry_s[ pos_y < 0.5 ]
e_errz_s1 = errz_s[ pos_z < 0.5 ]

e_errx_s2 = errx_s[ pos_x > 0.5 ]
e_erry_s2 = erry_s[ pos_y > 0.5 ]
e_errz_s2 = errz_s[ pos_z > 0.5 ]

e_errx_bh1 = errx_bh[ pos_x < 0.5 ]
e_erry_bh1 = erry_bh[ pos_y < 0.5 ]
e_errz_bh1 = errz_bh[ pos_z < 0.5 ]

e_errx_bh2 = errx_bh[ pos_x > 0.5 ]
e_erry_bh2 = erry_bh[ pos_y > 0.5 ]
e_errz_bh2 = errz_bh[ pos_z > 0.5 ]



figure(frameon=True)
subplot(311, title="Acceleration along X")
hist(e_errx_s1, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r', label="Sorted")
hist(e_errx_s2, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r', ls='dashed')
hist(e_errx_bh1, bins=bins, normed=1, histtype='step', rwidth=0.01, color='b', label="B-H")
hist(e_errx_bh2, bins=bins, normed=1, histtype='step', rwidth=0.01, color='b', ls='dashed')
legend(loc="upper right")
xlim(-0.05, 0.05)

subplot(312, title="Acceleration along Y")
hist(e_erry_s1, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r')
hist(e_erry_s2, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r', ls='dashed')
hist(e_erry_bh1, bins=bins, normed=1, histtype='step', rwidth=0.01, color='b', label="B-H")
hist(e_erry_bh2, bins=bins, normed=1, histtype='step', rwidth=0.01, color='b', ls='dashed')
xlim(-0.05, 0.05)

subplot(313, title="Acceleration along Z")
hist(e_errz_s1, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r')
hist(e_errz_s2, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r', ls='dashed')
hist(e_errz_bh1, bins=bins, normed=1, histtype='step', rwidth=0.01, color='b', label="B-H")
hist(e_errz_bh2, bins=bins, normed=1, histtype='step', rwidth=0.01, color='b', ls='dashed')
xlim(-0.05, 0.05)

savefig("histogram_all.png")
close()




print "E_error for sorted: ( x= %4.3e"%std(e_errx_s), "y= %4.3e"%std(e_erry_s), "z= %4.3e"%std(e_errz_s), ") Combined YZ = %4.3e"%(( std(e_erry_s) + std(e_errz_s) )/2.)

    #print std(e_errx_s), std(e_erry_s), std(e_errz_s)


    #print "Min   for sorted: ( x= %4.3e"%min(errx_s), "y= %4.3e"%min(erry_s), "z= %4.3e"%min(errz_s), ") Combined YZ = %4.3e"%(min( min(erry_s), min(errz_s) ))
    #print "Max   for sorted: ( x= %4.3e"%max(errx_s), "y= %4.3e"%max(erry_s), "z= %4.3e"%max(errz_s), ") Combined YZ = %4.3e"%(max( max(erry_s), max(errz_s) ))



