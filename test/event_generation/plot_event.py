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

def main(args):

    R=5

    data = np.loadtxt("event0_freq.dat")

    tmin=0.05
    tmax=0.1
    wmin=24.6
    wmax=26.2

    data = data[data[:,1]!=0]
    fig, ax = plt.subplots()
    #ax.plot(data[:,0], data[:,1]/(2*np.pi*1e9), ls='None', marker=',')
    ax.plot(data[:,0], data[:,1]/(2*np.pi), ls='None', marker=',')
    #ax.hist2d(data[:,0], data[:,1]/(2*np.pi),bins=10000, range=[[0,tmax],[0,wmax]])
    ax.set_xlabel('t[s]')
    ax.set_ylabel('f[GHz]')
    ax.set_xlim([tmin,tmax])
    ax.set_ylim([wmin,wmax])
    plt.show()

    data = np.loadtxt("event0_pos.dat")

    fig, ax = plt.subplots()
    #ax.plot(data[:,0], data[:,1]/(2*np.pi*1e9), ls='None', marker=',')
    ax.hist2d(data[:,0], data[:,1],bins=100,range=[[-R,R],[-R,R]])
    #ax.scatter(data[:,0], data[:,1], marker='.')
    ax.set_aspect('equal')
    ax.set_xlabel('x[cm]')
    ax.set_ylabel('y[cm]')
    ax.set_xlim([-R,R])
    ax.set_ylim([-R,R])
    plt.show()

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
