#!/usr/bin/env python3
# coding:utf-8
import matplotlib
import matplotlib.pyplot as plt # v1.5.3
import matplotlib.ticker as ticker
import matplotlib.animation as animation
import numpy as np
import sys
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
  
if __name__ == "__main__":

    data = np.genfromtxt(sys.argv[1])
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    u = data[:,3] + 1j*data[:,4]

    # cs = np.real(u)
    # vmin = min(cs)
    # vmax = max(cs)

    cs = np.abs(u)
    vmin = 0.0
    vmax = max(cs)

    cm = plt.get_cmap("jet")
    cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x, y, z, c=scalarMap.to_rgba(cs))
    scalarMap.set_array(cs)
    
    
    cb = fig.colorbar(scalarMap)
    
    plt.show()
    
    # figを閉じる
    fig.clf()
    plt.close(fig) # これがないと完全にfig,axがresetされない
