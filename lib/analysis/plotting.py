#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  plot_result.py
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
import argparse
import os
import seaborn as sns

sns.set()
        
class Plot_Wrapper:

    def __init__(self):

        #self.fig, self.ax = plt.subplots()
        self.initialized = False

    def add(self, f, *args):

        if not self.initialized:
            self.fig, self.ax = plt.subplots()
            
        f(self.fig, self.ax, *args)

    def finish(self, name=''):

        if name=='':
            display(self.fig)
        else:
            self.fig.savefig(name)
        
        plt.close(self.fig)
        self.initialized = False
        
    def save(self, name):
        
        self.fig.savefig(name)

    @classmethod
    def plot(cls, f, *args, name=''):

        #just a shorthand for quick plots

        tmp = cls()
        tmp.add(f, *args)
        tmp.finish(name)
            

def beamforming(fig, ax, R, data, cbar_label='', ax_unit='cm'):
    
    im_masked = np.ma.masked_where(data==0,data)
    im=ax.imshow(im_masked,extent=(-R,R,-R,R),origin='lower', zorder=2)

    cbar = fig.colorbar(im)
    ax.set_aspect('equal')
    ax.set_xlim(-(R+0.5),R+0.5)
    ax.set_ylim(-(R+0.5),R+0.5)
    
    xlabel='x[' + ax_unit + ']'
    ylabel='y[' + ax_unit + ']'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    cbar.ax.set_ylabel(cbar_label)


def main(args):
	
    binary = "../../bin/reconstruction "
    R=5
    r=3.5
    w0=26.0016
    phi = 0
    snr=0.5
    seed=-1
    n_samples=1000
    grid_size=100
    N = 30

    data = run_reconstruction(binary, grid_size, n_samples, snr, seed, r, phi, w0, N)
    
    plot_result(R, data[438], "beamforming_rec_test.pdf")

    data = np.loadtxt("beamforming_rec_ref.dat")
    plot_result(R, data, "beamforming_rec_ref.pdf")

    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
