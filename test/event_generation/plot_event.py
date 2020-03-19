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

    data = np.loadtxt("event0.dat")

    print(data.shape)
    data = data[data[:,1]!=0]
    fig, ax = plt.subplots()
    #ax.plot(data[:,0], data[:,1]/(2*np.pi*1e9), ls='None', marker=',')
    ax.plot(data[:,0], data[:,1]/(2*np.pi), ls='None', marker=',')
    ax.set_xlabel('t[s]')
    ax.set_ylabel('f[GHz]')
    #ax.set_xlim([0,1])
    #ax.set_ylim([0,1.6])
    plt.show()

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
