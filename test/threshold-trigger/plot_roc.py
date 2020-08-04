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

def plot_curve(ax, data, label):

    ind = np.argsort(data[:,0])
    x = data[ind,0]
    y = data[ind,2]
    auc = np.trapz(y=y, x=x)
    ax.errorbar(data[:,0], data[:,2], yerr=data[:,3], xerr=data[:,1], 
                                label=label + " AUC={:.2f}".format(auc), ls='None')

def main(args):

    sns.set()
    plt.rcParams.update({'lines.markeredgewidth': 1})
    
    n_samples = 1000
    grid_size_vals = [25, 50, 100]
    snr_vals = [0.1, 0.5, 1]
    n_packets = 5000
    
    t = int( time.time() * 1000.0 )
    seed = ((t & 0xff000000) >> 24) +\
             ((t & 0x00ff0000) >>  8) +\
             ((t & 0x0000ff00) <<  8) +\
             ((t & 0x000000ff) << 24)
    print(seed)
    
    seed = 123456

    i=0
    for snr in snr_vals:
        fig, ax = plt.subplots()
        fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.68,
					wspace=0.01, hspace=0.1)
        for grid_size in grid_size_vals:
            os.system("../../bin/threshold-trigger " + str(grid_size) + " "\
                        + str(n_samples) + " " + str(snr) + " "\
                        + str(n_packets) + " " + str(seed) +\
                         " >result_" + str(i) + ".out")

            data = np.loadtxt("result_"+str(i)+".out")
            #plt.plot(data[:,0], data[:,1])#, ls='None', marker='o')

            plot_curve(ax, data, label="g="+str(grid_size))

            os.system("rm result_"+str(i)+".out")
			
            i+=1
			
        ax.set_title("snr="+str(snr))
        ax.set_ylabel("sensitivity")
        ax.set_xlabel("1-specificity")
        ax.set_xlim(-0.1,1.1)
        ax.set_ylim(-0.1,1.1)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig("roc_"+str(int(10*snr))+".pdf")
		
		
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
