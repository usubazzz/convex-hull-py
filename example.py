#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np
import numpy.random as rd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

import convexhull

def plot_ridge_2D(facets, points, ax, c):
    ims = []
    for f in facets:
        line_x = [points[pid[0]][0] for pid in f.ridge]
        line_y = [points[pid[0]][1] for pid in f.ridge]
        ims += ax.plot(line_x, line_y, "-", color=c)
    return ims

def plot_ridge_3D(facets, points, ax, c):
    ims = []
    for f in facets:
        for r in f.ridge:
            line_x = [points[pid][0] for pid in r]
            line_y = [points[pid][1] for pid in r]
            line_z = [points[pid][2] for pid in r]
            ims += ax.plot(line_x, line_y, line_z, "-", color=c)
    return ims

def plot_normal_3D(facets, points, ax, c):
    ims = []
    for f in facets:
        centroid = np.array([0.0, 0.0, 0.0])
        for pid in f.points_id:
            centroid += np.array(points[pid])
        offset = centroid / 3.0
        ims += ax.plot([offset[0], offset[0] + f.normal[0]], [offset[1], offset[1] + f.normal[1]], [offset[2], offset[2] + f.normal[2]], "-", color=c)

    return ims


# ランダム範囲 min <= rnd < max
rnd_min = -10
rnd_max = 10
seed = rd.randint(0, 300)
rd.seed(seed)
print("seed: {}".format(seed))

N = 30  # データ数
D = 3   # データ次元数
points = [[(rnd_max-rnd_min)*rd.random() + rnd_min for i in range(D)] for j in range(N)] # 計算用の頂点群

body = convexhull.QuickHull(points, D)

body.first_hull()

fig = plt.figure(figsize=plt.figaspect(1.0))
if body.dim == 2:
    ax = fig.add_subplot(1,1,1)
if body.dim == 3:
    ax = Axes3D(fig)

ims = []

### 頂点表示
if body.dim == 2:
    for v in points:
        ax.scatter(v[0], v[1], marker="o", c='r', s=10)
if body.dim == 3:
    for v in points:
        ax.scatter(v[0], v[1], v[2], marker="o", c='r', s=10)

### 初回ファセット
if body.dim == 2:
    ims1 = plot_ridge_2D(body.facets, body.points, ax, 'k')
    ims2 = plot_ridge_2D(body.stacks, body.points, ax, 'b')
    imc = ax.scatter(body.centroid[0], body.centroid[1], marker="o", c='g', s=60)
    ims.append(ims1+ims2+[imc])
if body.dim ==3:
    ims1 = plot_ridge_3D(body.facets, body.points, ax, 'k')
    ims2 = plot_ridge_3D(body.stacks, body.points, ax, 'b')
    ims3 = plot_normal_3D(body.stacks, body.points, ax, 'g')
    ims4 = plot_normal_3D(body.facets, body.points, ax, 'g')
    imc = ax.scatter(body.centroid[0], body.centroid[1], body.centroid[2], marker="o", c='g', s=60)
    ims.append(ims1+ims2+ims3+ims4+[imc])

### メインループ
while body.stacks != []:
    _facet = body.stacks[-1]

    body.make_hull_step(_facet)

    if body.dim == 2:
        ims1 = plot_ridge_2D(body.facets, body.points, ax, 'k')
        ims2 = plot_ridge_2D(body.stacks, body.points, ax, 'b')
        imc = ax.scatter(body.centroid[0], body.centroid[1], marker="o", c='g', s=60)
        ims.append(ims1+ims2+[imc])
    if body.dim == 3:
        ims1 = plot_ridge_3D(body.facets, body.points, ax, 'k')
        ims2 = plot_ridge_3D(body.stacks, body.points, ax, 'b')
        ims3 = plot_normal_3D(body.stacks, body.points, ax, 'g')
        ims4 = plot_normal_3D(body.facets, body.points, ax, 'g')
        imc = ax.scatter(body.centroid[0], body.centroid[1], body.centroid[2], marker="o", c='g', s=60)
        ims.append(ims1+ims2+ims3+ims4+[imc])

print("seed: {}".format(seed))



plt.xlim([rnd_min, rnd_max])
plt.ylim([rnd_min, rnd_max])

ani = animation.ArtistAnimation(fig, ims, interval=1000)

# ani.save("sample.gif", writer="imagemagick")

plt.show()