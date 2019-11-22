#!/usr/bin/env python3
# coding:utf-8
from stl import mesh
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np

def show(data):
    figure = plt.figure()
    axes = mplot3d.Axes3D(figure)

    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(data.vectors))
    scale = data.points.flatten(-1)
    axes.auto_scale_xyz(scale, scale, scale)
    plt.show()

  
if __name__ == "__main__":
    data = mesh.Mesh.from_file('sphere.stl')
    
    show(data)
    
