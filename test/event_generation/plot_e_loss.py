#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  plot_e_loss.py
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
from scatter_dis import plot_loss_dis

def main(args):

    data = np.loadtxt('e_loss.dat')

    data_r=data[data<50]

    fig, ax = plt.subplots()

    plot_loss_dis(ax)
    ax.hist(data_r, density=True, bins=100)
    fig.savefig("e_loss.pdf")
    plt.close(fig)

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
