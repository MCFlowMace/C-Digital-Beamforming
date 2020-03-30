#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  scatter_dis.py
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

c=14.0921
b=2.5*c

A1 = 0.204
w1 = 1.85
e1 = 12.6
A2 = 0.0556
w2 = 12.5
e2 = 14.3

import numpy as np
import matplotlib.pyplot as plt

def f1(x,A,e,w):
    return A*np.exp(-2*(x-e)**2/w**2)

def f2(x,A,e,w):
    return A*w**2/(w**2 + 4*(x - e)**2)

def f(x):

    #~ if(x<c):
        #~ return f1(x,A1,e1,w1)
    #~ else:
        #~ return f2(x, A2, e2, w2)
    return np.where(x<c, f1(x,A1,e1,w1), f2(x, A2, e2, w2))

def Finv(x, t, s):
    return t + s*np.tan(0.5*np.pi*(2*x-1))

def bw(n):
    u=np.random.uniform(0.5,1.0,n)
    return Finv(u, e2, w2/2)

def triangle(n):
    u=np.random.rand(n)

    vals = b-np.sqrt((1-u)*(b-c)**2)
    return vals

def replace(data, c):
    l = int(data.size*c)
    data_m = data
    data_m[:l]=bw(l)
    return data_m

def main(args):

    x=np.linspace(0,50,1000)
    data_f=f(x)

    data = np.random.normal(e1,w1/2,100000)
    data_m=replace(data,0.5)
    data_r=data_m[(data_m>0)&(data_m<50)]
    plt.hist(data_r, density=True,bins=100)
    plt.plot(x,data_f)
    plt.show()

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
