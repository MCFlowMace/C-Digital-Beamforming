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
import io
import subprocess
from pathlib import Path


def plot_result(R, data, name, save=True):

    fig, ax = plt.subplots()

    im_masked = np.ma.masked_where(data==0,data)
    im = ax.imshow(im_masked.T, extent=(-R,R,-R,R), origin='lower')
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
    
    if save:
        plt.savefig(name)
    else:
        plt.show()
        
    plt.close(fig)


def run_reconstruction(binary):
    result = subprocess.run(
        [str(binary)],
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    )

    # Beamformer prints this diagnostic before the matrix.
    output = "\n".join(
        line for line in result.stdout.splitlines()
        if line.strip() != "running run_next"
    )

    data = np.loadtxt(io.StringIO(output))
    data[~np.isfinite(data)] = 0
    data[data == -1] = 0

    return data


def main(args):
    binary = Path(__file__).resolve().parents[2] / "build/bin/test_rec"
    R = 5

    data = run_reconstruction(binary)
    plot_result(R, data, "beamforming_rec_test.pdf")

    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
