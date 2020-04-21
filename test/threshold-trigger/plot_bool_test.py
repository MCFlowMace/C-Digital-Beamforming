#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  plot_bool_test.py
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

    data = np.loadtxt('output.dat')

    fig, ax1 = plt.subplots()

    print(data[data[:,0]!=0,0])

    ax1.plot(data[:,0], label='frequency', ls='None', c='r', marker='x', ms=0.1)

    ax2 = ax1.twinx()
    ax2.plot(data[:,1], label='label', ls='None', c='b', marker='+', ms=0.1)
    ax2.set_ylim([0,10])

    ax1.legend(loc='best')
    ax2.legend(loc='best')
    plt.show()

    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
