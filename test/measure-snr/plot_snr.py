#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  plot_snr.py
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

def main(args):

    parser = argparse.ArgumentParser(description='Plots the snr result')
    parser.add_argument('input', metavar='path', type=str,
                       help='path for the input')
    parser.add_argument('mode', metavar='mode', type=str,
                        help='0 plots against samples 1 against antennas')
    args = parser.parse_args()
    inFile = args.input
    mode = args.mode

    print(mode)

    if mode=='True':
        print("samples")
        label='samples'
        ind_x=1
    else:
        label='N'
        ind_x=0

    data = np.loadtxt(inFile)

    fig, ax = plt.subplots()

    ax.errorbar(data[:,ind_x], data[:,2], yerr=data[:,5],ls='None', marker='.', capsize=3.0)
    ax.set_xlabel(label)
    ax.set_ylabel('SNR')

    fig.savefig(inFile+'.pdf')
    plt.close(fig)


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
