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
from scipy import stats

dist_cutoff_ratio=1.2

print "Plotting..."

axis = [
    1.0,  1.0,  1.0, 
    1.0,  1.0,  0.0, 
    1.0,  1.0, -1.0, 
    1.0,  0.0,  1.0, 
    1.0,  0.0,  0.0, 
    1.0,  0.0, -1.0, 
    1.0, -1.0,  1.0, 
    1.0, -1.0,  0.0, 
    1.0, -1.0, -1.0, 
    0.0,  1.0,  1.0, 
    0.0,  1.0,  0.0, 
    0.0,  1.0, -1.0, 
    0.0,  0.0,  1.0,
    -1.0, -1.0, -1.0, 
    -1.0, -1.0,  0.0, 
    -1.0, -1.0,  1.0, 
    -1.0,  0.0, -1.0, 
    -1.0,  0.0,  0.0, 
    -1.0,  0.0,  1.0, 
    -1.0,  1.0, -1.0, 
    -1.0,  1.0,  0.0, 
    -1.0,  1.0,  1.0, 
     0.0, -1.0, -1.0, 
     0.0, -1.0,  0.0, 
     0.0, -1.0,  1.0, 
     0.0,  0.0, -1.0  
  ]


#names = ["side", "edge", "corner"]

#for orientation in range( 26 ):
for jjj in range(3):
    if jjj == 0:
        orientation = 0
    if jjj == 1:
        orientation = 1
    if jjj == 2:
        orientation = 4

         
    # Read Quickshed accelerations
    data=loadtxt( "interaction_dump_%d.dat"%orientation )
    id = data[:,0]
    pos = data[:,2]
    pos -= mean(pos)

    x = data[:,3]
    y = data[:,4]
    z = data[:,5]

    dist = data[:,12]
    
    accx_u=data[:,6]
    accy_u=data[:,7]
    accz_u=data[:,8]
    
    accx_s=data[:,9]
    accy_s=data[:,10]
    accz_s=data[:,11]
    



    # Build error ------------------------------------------------
    
    inv_acc_tot = 1.0 / sqrt(accx_u**2 + accy_u**2 + accz_u**2) 
    errx_s = (accx_s - accx_u ) * inv_acc_tot
    erry_s = (accy_s - accy_u ) * inv_acc_tot
    errz_s = (accz_s - accz_u ) * inv_acc_tot
    
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
    

#     sample_pos = pos[dist<0.2]
#     sample_x = e_errx_s[dist<0.2]
#     sample_y = e_erry_s[dist<0.2]
#     sample_z = e_errz_s[dist<0.2]
            

#     numBins = 100
#     binEdges = linspace(-1.5-1.5/numBins, 1.5+1.5/numBins, numBins+2)
#     bins = linspace(-1.5, 1.5, numBins+1)

#     corr_x, a, b = stats.binned_statistic(sample_pos, sample_x, statistic='mean', bins=binEdges)
#     corr_y, a, b = stats.binned_statistic(sample_pos, sample_y, statistic='mean', bins=binEdges)
#     corr_z, a, b = stats.binned_statistic(sample_pos, sample_z, statistic='mean', bins=binEdges)

#     a,b, sample_bin = stats.binned_statistic(pos, pos, statistic='mean', bins=binEdges)
#     sample_bin -= 2

    
# #    for j in range(size(pos)):
# #        e_errx_s /= corr_x[sample_bin[j]]
# #        e_erry_s /= corr_y[sample_bin[j]]
# #        e_errz_s /= corr_z[sample_bin[j]]
            
    # Plot -------------------------------------------------------
    figure(frameon=True)

    
    subplot(311, title="Acceleration along X")
    scatter(pos[dist<1.], e_errx_s[dist<1.] ,c=dist[dist<1.], marker='o', s=1, linewidth=0, cmap='jet')
    plot([-dist_cutoff_ratio/2., -dist_cutoff_ratio/2.],  [-0.01, 0.01], 'k--')
    plot([dist_cutoff_ratio/2., dist_cutoff_ratio/2.],  [-0.01, 0.01], 'k--')
    text( 0., 0.005, "axis=( %d %d %d )"%(axis[orientation*3 + 0], axis[orientation*3 + 1], axis[orientation*3 + 2]) , ha='center', backgroundcolor='w', fontsize=14)
    xlim(-1.2*max(abs(pos)), 1.2*max(abs(pos)))
    ylim(-0.03, 0.03)
    grid()
    
    subplot(312, title="Acceleration along Y")
    scatter(pos[dist<1.], e_erry_s[dist<1.] , c=dist[dist<1.], marker='o', s=1, linewidth=0, cmap='jet')
    plot([-dist_cutoff_ratio/2., -dist_cutoff_ratio/2.],  [-0.03, 0.03], 'k--')
    plot([dist_cutoff_ratio/2., dist_cutoff_ratio/2.],  [-0.03, 0.03], 'k--')
    text( 0., 0.005, "axis=( %d %d %d )"%(axis[orientation*3 + 0], axis[orientation*3 + 1], axis[orientation*3 + 2]) , ha='center', backgroundcolor='w', fontsize=14)
    xlim(-1.2*max(abs(pos)), 1.2*max(abs(pos)))
    ylim(-0.03, 0.03)  
    grid()
    
    subplot(313, title="Acceleration along Z")
    scatter(pos[dist<1.], e_errz_s[dist<1.] , c=dist[dist<1.], label="Sorted", marker='o', s=1, linewidth=0, cmap='jet')
    plot([-dist_cutoff_ratio/2., -dist_cutoff_ratio/2.],  [-0.03, 0.03], 'k--')
    plot([dist_cutoff_ratio/2., dist_cutoff_ratio/2.],  [-0.03, 0.03], 'k--')
    text( 0., 0.005, "axis=( %d %d %d )"%(axis[orientation*3 + 0], axis[orientation*3 + 1], axis[orientation*3 + 2]) , ha='center', backgroundcolor='w', fontsize=14)
    #legend(loc="upper right")
    xlim(-1.2*max(abs(pos)), 1.2*max(abs(pos)))
    ylim(-0.03, 0.03)
    grid()

    savefig("quadrupole_accelerations_relative_%d.png"%orientation)
    close()



    figure(frameon=True)
    
    subplot(311, title="Acceleration along X")
    scatter(dist, e_errx_s )
    text( 0., 0.005, "axis=( %d %d %d )"%(axis[orientation*3 + 0], axis[orientation*3 + 1], axis[orientation*3 + 2]) , ha='center', backgroundcolor='w', fontsize=14)
    #xlim(-1.2*max(abs(pos)), 1.2*max(abs(pos)))
    #ylim(-0.01, 0.01)
    grid()
    
    subplot(312, title="Acceleration along Y")
    scatter(dist, e_erry_s )
    text( 0., 0.005, "axis=( %d %d %d )"%(axis[orientation*3 + 0], axis[orientation*3 + 1], axis[orientation*3 + 2]) , ha='center', backgroundcolor='w', fontsize=14)
    #xlim(-1.2*max(abs(pos)), 1.2*max(abs(pos)))
    #ylim(-0.01, 0.01)  
    grid()
    
    subplot(313, title="Acceleration along Z")
    scatter(dist, e_errz_s , label="Sorted")
    text( 0., 0.005, "axis=( %d %d %d )"%(axis[orientation*3 + 0], axis[orientation*3 + 1], axis[orientation*3 + 2]) , ha='center', backgroundcolor='w', fontsize=14)
    #legend(loc="upper right")
    #xlim(-1.2*max(abs(pos)), 1.2*max(abs(pos)))
    #ylim(-0.01, 0.01)
    grid()

    savefig("error_distance_%d.png"%orientation)
    close()


    
    # # Plot -------------------------------------------------------
    # figure(frameon=True)
    
    # subplot(311, title="Acceleration along X")
    # #plot(id[abs(errx_s) > 0.001], e_errx_s , 'ro')
    # plot(pos, accx_u , 'bx')
    # plot(pos, accx_s , 'ro')
    # text( 0., 0.1, "axis=( %d %d %d )"%(axis[orientation*3 + 0], axis[orientation*3 + 1], axis[orientation*3 + 2]) , ha='center', backgroundcolor='w', fontsize=14)
    # ylim(-700, 700)
    # grid()
    
    # subplot(312, title="Acceleration along Y")
    # #plot(id[abs(erry_s) > 0.001], e_erry_s , 'ro')
    # plot(pos, accy_u , 'bx')
    # plot(pos, accy_s , 'ro')
    # text( 0., 0.1, "axis=( %d %d %d )"%(axis[orientation*3 + 0], axis[orientation*3 + 1], axis[orientation*3 + 2]) , ha='center', backgroundcolor='w', fontsize=14)
    # ylim(-700, 700) 
    # grid()
    
    # subplot(313, title="Acceleration along Z")
    # #plot(id[abs(errz_s) > 0.001], e_errz_s , 'ro', label="Sorted")
    # plot(pos, accz_u , 'bx', label="Unsorted")
    # plot(pos, accz_s , 'ro', label="Sorted")

    # text( 0., 0.1, "axis=( %d %d %d )"%(axis[orientation*3 + 0], axis[orientation*3 + 1], axis[orientation*3 + 2]) , ha='center', backgroundcolor='w', fontsize=14)

    # legend(loc="upper right")
    
    # ylim(-700, 700)
    # grid()

    # savefig("accelerations_absolute_%d.png"%orientation)
    # close()
    
    
    
    # # Plot -------------------------------------------------------
    # bins = linspace(-3, 3, 10000)
    
    # e_errx_s = errx_s[(abs(errx_s) >= 0.000) & (abs(errx_s) < 1.)]
    # e_erry_s = erry_s[(abs(erry_s) >= 0.000) & (abs(erry_s) < 1.)]
    # e_errz_s = errz_s[(abs(errz_s) >= 0.000) & (abs(errz_s) < 1.)]
    
    
    
    figure(frameon=True)
    subplot(311, title="Acceleration along X")
    hist(e_errx_s, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r', label="Sorted")
    legend(loc="upper right")
    xlim(-0.12, 0.12)
    
    subplot(312, title="Acceleration along Y")
    hist(e_erry_s, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r')
    xlim(-0.12, 0.12)

    subplot(313, title="Acceleration along Z")
    hist(e_errz_s, bins=bins, normed=1, histtype='step', rwidth=0.01, color='r')
    xlim(-0.12, 0.12)
    
    savefig("histogram_%d.png"%orientation)
    close()




    print "E_error for sorted: ( x= %4.3e"%std(e_errx_s), "y= %4.3e"%std(e_erry_s), "z= %4.3e"%std(e_errz_s), ") Combined YZ = %4.3e"%(( std(e_erry_s) + std(e_errz_s) )/2.)

    #print std(e_errx_s), std(e_erry_s), std(e_errz_s)


    #print "Min   for sorted: ( x= %4.3e"%min(errx_s), "y= %4.3e"%min(erry_s), "z= %4.3e"%min(errz_s), ") Combined YZ = %4.3e"%(min( min(erry_s), min(errz_s) ))
    #print "Max   for sorted: ( x= %4.3e"%max(errx_s), "y= %4.3e"%max(erry_s), "z= %4.3e"%max(errz_s), ") Combined YZ = %4.3e"%(max( max(erry_s), max(errz_s) ))



 
