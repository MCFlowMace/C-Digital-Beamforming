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
from scipy.fft import rfft, rfftfreq, fft, fftfreq, irfft, fftshift, ifft, fft2
import os

def plot_result(R, data, name):

    fig, ax = plt.subplots()

    im_masked = np.ma.masked_where(data==-1,data)
    im=ax.imshow(im_masked,extent=(-R,R,-R,R),origin='lower')
    #im=ax.imshow(np.transpose(self.grid_phis[0]),extent=(-R,R,-R,R),origin='lower')
    #ax.plot(electron.x,electron.y,c='r',marker='o',ms=2)
    fig.colorbar(im)
    ax.set_aspect('equal')
    ax.set_xlim(-(R+0.5),R+0.5)
    ax.set_ylim(-(R+0.5),R+0.5)
    #ax.axhline(0, -R, R, c='b')
    xlabel='x[cm]'
    ylabel='y[cm]'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.savefig(name)
    plt.close(fig)

def run_reconstruction(binary, grid_size, n_samples, snr, 
							seed, r, phi, w0, N, weighted=1):
                                
    cmd = binary + str(grid_size) + " "\
                + str(n_samples) + " " + str(snr) + " "\
                + str(seed) + " " + str(weighted) + " "\
                + str(r) + " " + str(phi) + " " + str(w0) + " "\
                + str(N) + " >result.out"
    
    os.system(cmd)
                
    data = np.loadtxt("result.out")
    os.system("rm result.out")
    
    data = data.reshape((-1, grid_size, grid_size))
    
    return data

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
    
    print(data.shape)
    plot_result(R, data[438], "beamforming_rec_test.pdf")
    
    data[data==-1]=0
    
    ind = (np.isfinite(data))^True
    data[ind] = 0
    plot_result(R, data, "beamforming_rec_test.pdf")
    
    data_masked = np.ma.masked_where(data==-1,data)
    data_freq = fftshift(fft2(data_masked))
    
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    im = ax.imshow(np.transpose(np.abs(data_freq)), origin='lower')
    fig.colorbar(im)
    plt.savefig("beamforming_freq_test.pdf")
    plt.close(fig)


    data = np.loadtxt("beamforming_rec_ref.dat")
    plot_result(R, data, "beamforming_rec_ref.pdf")

    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
