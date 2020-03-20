#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  plot_event.py
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

def plot_freq(tmin, tmax, fmin, fmax, data, name):

    fig, ax = plt.subplots()
    #ax.plot(data[:,0], data[:,1]/(2*np.pi*1e9), ls='None', marker=',')
    ax.plot(data[:,0], data[:,1]/(2*np.pi), ls='None', marker=',')
    #ax.hist2d(data[:,0], data[:,1]/(2*np.pi),bins=10000, range=[[0,tmax],[0,wmax]])
    ax.set_xlabel('t[s]')
    ax.set_ylabel('f[GHz]')
    ax.set_xlim([tmin,tmax])
    ax.set_ylim([fmin,fmax])
    plt.savefig(name, dpi=600)
    plt.close(fig)

def plot_pos(R, data, name):

    fig, ax = plt.subplots()
    #ax.plot(data[:,0], data[:,1]/(2*np.pi*1e9), ls='None', marker=',')
    ax.hist2d(data[:,0], data[:,1],bins=100,range=[[-R,R],[-R,R]])
    #ax.scatter(data[:,0], data[:,1], marker='.')
    ax.set_aspect('equal')
    ax.set_xlabel('x[cm]')
    ax.set_ylabel('y[cm]')
    ax.set_xlim([-R,R])
    ax.set_ylim([-R,R])
    plt.savefig(name, dpi=300)
    plt.close(fig)

def main(args):

    parser = argparse.ArgumentParser(description='Plots the event')
    parser.add_argument('input1', metavar='path1', type=str,
                       help='path for the frequency input')
    parser.add_argument('tmin', metavar='tmin', type=float,
                        help='lower bound of time window')
    parser.add_argument('tmax', metavar='tmax', type=float,
                        help='upper bound of time window')
    parser.add_argument('fmin', metavar='fmin', type=float,
                        help='lower bound of frequency window')
    parser.add_argument('fmax', metavar='fmax', type=float,
                        help='upper bound of frequency window')

    parser.add_argument('input2', metavar='path2', type=str,
                        help='path for the position input')
    parser.add_argument('R', metavar='radius', type=float,
                        help='Radius for the position plot')


    args = parser.parse_args()

    inFile1 = args.input1
    tmin = args.tmin
    tmax = args.tmax
    fmin = args.fmin
    fmax = args.fmax

    inFile2 = args.input2
    R=args.R

    data = np.loadtxt(inFile1)
    data = data[data[:,1]!=0]
    plot_freq(tmin, tmax, fmin, fmax, data, "event0_spec.pdf")

    data = np.loadtxt(inFile2)
    plot_pos(R, data, "event0_pos.png")


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
