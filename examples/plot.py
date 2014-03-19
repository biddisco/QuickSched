#!/usr/bin/env python2                                                                                                                     
# -*- coding: utf-8 -*-                                                                                                                    
import matplotlib
matplotlib.use("Agg")
from pylab import *
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

data=loadtxt("particle_dump.dat")
id = data[:,0]
accx_e=data[:,1]
accy_e=data[:,2]
accz_e=data[:,3]

accx_bh=data[:,4]
accy_bh=data[:,5]
accz_bh=data[:,6]

accx_new=data[:,7]
accy_new=data[:,8]
accz_new=data[:,9]



figure(frameon=True)

subplot(311)
plot(id, (accx_bh - accx_e )/abs(accx_e) + 1, 'rx')
plot(id, (accx_new - accx_e )/abs(accx_e)  + 1 , 'b.')

ylim(0.8, 1.2)
grid()

subplot(312)
plot(id, (accy_bh - accy_e )/abs(accy_e ) + 1, 'rx')
plot(id, (accy_new - accy_e )/abs(accy_e ) + 1, 'b.')
ylim(0.8, 1.2)
grid()

subplot(313)
plot(id, (accz_bh - accz_e )/abs(accz_e ) + 1, 'rx')
plot(id, (accz_new - accz_e )/abs(accz_e)  + 1, 'b.')
ylim(0.8, 1.2)
grid()

savefig("accelerations.png")

