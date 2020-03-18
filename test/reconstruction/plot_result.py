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

def plot_result(R, data, name):

    fig, ax = plt.subplots()

    im_masked = np.ma.masked_where(data==-1,data)
    im=ax.imshow(np.transpose(im_masked),extent=(-R,R,-R,R),origin='lower')
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


def main(args):

    parser = argparse.ArgumentParser(description='Plots the beamforming result')
    parser.add_argument('input', metavar='path', type=str,
                       help='path for the input')
    parser.add_argument('R', metavar='radius', type=float,
                        help='Radius for the plot')
    args = parser.parse_args()
    inFile = args.input
    R=args.R

    data = np.loadtxt(inFile)
    plot_result(R, data, "beamforming_rec_test.pdf")

    data = np.loadtxt("beamforming_rec_ref.dat")
    plot_result(R, data, "beamforming_rec_ref.pdf")

    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
