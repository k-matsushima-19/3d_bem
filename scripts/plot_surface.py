from stl import mesh
import matplotlib
import mpl_toolkits.axes_grid1
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.cm as cmx
import numpy as np
import sys

# 参考: https://scipython.com/book/chapter-8-scipy/examples/visualizing-the-spherical-harmonics/
def plot_surface(fig, ax, file_stl, file_sol, mode="real"):

    data = mesh.Mesh.from_file(file_stl)

    # polygonを生成
    polys = mplot3d.art3d.Poly3DCollection(data.vectors, lw=0)
    scale = data.points.flatten(-1)

    # 各要素の中点
    centers = []
    for verts in data.vectors:
        centers.append(np.sum(verts,axis=0))
    centers = np.array(centers)
    
    # 各要素の色
    fcolors = centers[:,0]
    data = np.genfromtxt(file_sol)
    u = data[:,3] + 1j*data[:,4]

    if mode == "real":
        fcolors = np.real(u)
    elif mode == "imag":
        fcolors = np.imag(u)
    elif mode == "abs":
        fcolors = np.abs(u)

    # [-1,1]に正規化
    fmax, fmin = fcolors.max(), fcolors.min()
    # fcolors = (fcolors - fmin)/(fmax - fmin)
    fcolors_normalized = (fcolors - fmin)/(fmax - fmin)
    
    # 色を塗る
    # polys.set_facecolor(cm.seismic(fcolors))
    polys.set_facecolor(cm.jet(fcolors_normalized))
    
    image = ax.add_collection3d(polys)

    # aspect比を揃える
    ax.auto_scale_xyz(scale, scale, scale)

    # color bar
    vmin = min(fcolors)
    vmax = max(fcolors)
    cmap = plt.get_cmap("jet")
    cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
    scalarMap.set_array(fcolors)

    # divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
    # cax = divider.append_axes('right', '5%', pad='3%')
    cb = fig.colorbar(scalarMap, shrink=0.8)

if __name__ == "__main__":

    iwindow = 0
    for mode in ("real", "imag"):
        iwindow += 1
        fig = plt.figure(iwindow,figsize=(9,10))
        ax = mplot3d.Axes3D(fig)

        # plot_surface(fig, ax, "sphere.stl", "bndry.dat")

        plot_surface(fig, ax, sys.argv[1], sys.argv[2], mode=mode)

        ax.tick_params(pad=10)
        ax.tick_params(labelsize=15)
        ax.set_xlabel(r"$x_1$", fontsize=20, labelpad=15)
        ax.set_ylabel(r"$x_2$", fontsize=20, labelpad=15)
        ax.set_zlabel(r"$x_3$", fontsize=20, labelpad=15)

        
    
    plt.show()    


