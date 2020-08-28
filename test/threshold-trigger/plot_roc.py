#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  plot_roc.py
#
#  Copyright 2020 Florian Thomas <>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import time
from statsmodels.stats.proportion import proportion_confint

def get_roc(data):
	
	print("actual true: ", data[:,3])
	print("actual false: ", data[:,1])
	
	data_out = np.empty(shape=[data.shape[0], data.shape[1]+2])
	data_out[:,0] = data[:,0]/data[:,1]
	data_out[:,3] = data[:,2]/data[:,3]
	
	#x-error
	errors = proportion_confint(data[:,0], data[:,1])
	data_out[:,1] = data_out[:,0]-errors[0]
	data_out[:,2] = errors[1] - data_out[:,0]
	
	#y-error
	errors = proportion_confint(data[:,2], data[:,3])
	data_out[:,4] = data_out[:,3] - errors[0]
	data_out[:,5] = errors[1] - data_out[:,3]
	
	return data_out
    
def get_auc(data):
    
    ind = np.argsort(data[:,0])
    x = data[ind,0]
    y = data[ind,3]
    auc = np.trapz(y=y, x=x)
    
    return auc

def plot_curve(ax, data, label):
    
    auc = get_auc(data)

    ax.errorbar(data[:,0], data[:,3], xerr=[data[:,1], data[:,2]], 
								yerr=[data[:,4], data[:,5]], 
                                label=label + " AUC={:.2f}".format(auc), ls='None')
                                
    
def create_roc_plot(fig, ax, snr):

    ax.set_title("snr="+str(snr))
    ax.set_ylabel("sensitivity")
    ax.set_xlabel("1-specificity")
    ax.set_xlim(-0.1,1.1)
    ax.set_ylim(-0.1,1.1)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig("roc_"+str(int(10*snr))+".pdf")
        
def plot_curve_old(ax, data, label):

    ind = np.argsort(data[:,0])
    x = data[ind,0]
    y = data[ind,2]
    auc = np.trapz(y=y, x=x)
    ax.errorbar(data[:,0], data[:,2], xerr=data[:,1], yerr=data[:,3], 
                                label=label + " AUC={:.2f}".format(auc), ls='None')                               
                                
def run_threshold_trigger(binary, grid_size, n_samples, snr, 
                            n_packets, seed, r, phi, w0):
                                
    os.system(binary + str(grid_size) + " "\
                + str(n_samples) + " " + str(snr) + " "\
                + str(n_packets) + " " + str(seed) + " "\
                + str(r) + " " + str(phi) + " " + str(w0) + " "\
                + " >result.out")
                
    data = np.loadtxt("result.out")
    os.system("rm result.out")
    return data
    
def get_rand_seed():
    
    t = int( time.time() * 1000.0 )
    seed = ((t & 0xff000000) >> 24) +\
             ((t & 0x00ff0000) >>  8) +\
             ((t & 0x0000ff00) <<  8) +\
             ((t & 0x000000ff) << 24)
             
    return seed
    
def plot_x_AUC(x, AUC, snr_vals, name, xlabel, grid_boundaries):
    
    fig, ax = plt.subplots()
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.8,
                    wspace=0.01, hspace=0.1)
    
    for i in range(len(snr_vals)): 
        ax.plot(x, AUC[i], marker='^', label='snr='+str(snr_vals[i]))
    
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax)
    ax.vlines(grid_boundaries[grid_boundaries>=0],ymin=ymin, ymax=ymax, 
                                                ls='dashed', alpha=0.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("AUC")
    ax.set_title("grid-size="+str(grid_boundaries.shape[0]-1))
    
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.savefig(name+".pdf")
    plt.close(fig)

def main(args):

    sns.set()
    plt.rcParams.update({'errorbar.capsize': 1.5})
   
    R=5.0
    
    binary = "../../bin/threshold-trigger "
    n_samples = 1000
    grid_size_vals = [25, 50, 100]
    snr_vals = [0.1]
    n_packets = 1000
    r_vals = np.linspace(0, R, 100)
    phi = 0.0
    w0 = 26e9
    
    seed = get_rand_seed()
    print(seed)
    
    #seed = 123456
    
  
    grid_size = grid_size_vals[0]
    
    grid_boundaries = np.linspace(-R,R,grid_size+1)
    
    """
    
    auc_vals = np.empty(shape=[len(snr_vals),r_vals.shape[0]])
    
    for i in range(len(snr_vals)):
        for j, r in enumerate(r_vals):
        
            data = run_threshold_trigger(binary, grid_size_vals[0], n_samples, snr_vals[i], 
                                            n_packets, seed, r, phi, w0)
            
            data = get_roc(data)

            auc_vals[i,j] = get_auc(data) #plot_curve(ax, data, label="g="+str(grid_size_vals[0]))
    
    plot_x_AUC(r_vals, auc_vals, snr_vals, 'R_AUC', 'r[cm]', grid_boundaries)
        
   """

    for snr in snr_vals:
        fig, ax = plt.subplots()
        fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.68,
                    wspace=0.01, hspace=0.1)
        for grid_size in grid_size_vals:
            
            data = run_threshold_trigger(binary, grid_size, n_samples, snr, 
                            n_packets, seed, -1, phi, w0)
            
            data = get_roc(data)

            auc = plot_curve(ax, data, label="g="+str(grid_size))
            
        create_roc_plot(fig, ax, snr)
    
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
