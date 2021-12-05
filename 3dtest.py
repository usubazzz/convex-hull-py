import re
import numpy as np
from numpy.core.numeric import cross
from numpy.core.overrides import set_module
import numpy.random as rd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sympy

class Facet():
    def __init__(self, _id_list, _vertices):
        self.dim = len(_id_list)
        self.vertices = _vertices
        self.pints_id = _id_list
        self.edge = self.set_edge()
        self.normal = self.calc_normal()
        self.centroid = self.calc_centroid()

    def calc_normal(self):
        A = np.array([self.vertices[i] for i in sorted(self.pints_id)])
        x = np.dot(np.linalg.inv(A), np.ones(self.dim))         # pinv である必要はない？
        x = x / np.linalg.norm(x)
 
        p0 = self.vertices[self.pints_id[0]]
        p1 = self.vertices[self.pints_id[1]]
        # p2 = self.vertices[self.pints_id[2]]
        # p3 = self.vertices[self.pints_id[3]]

        u1 = np.array(p1) - np.array(p0)
        # u2 = np.array(p2) - np.array(p1)
        # u3 = np.array(p3) - np.array(p2)

        # ベクトルから求める方法
        #x = [-u1[1], u1[0]]

        print("normal(inv): {} dot: {}".format(x, np.round(np.dot(u1, x), 3)))
        # print("normal(inv): {} dot: {}".format(x, np.round(np.dot(u2, x), 3)))
        # print("normal(inv): {} dot: {}".format(x, np.round(np.dot(u3, x), 3)))

        return x

    def set_edge(self):
        point_list = sorted(self.pints_id)
        point_size = len(point_list)
        edge = []
        for i in range(point_size-1):
            edge.append([point_list[i], point_list[i+1]])
        edge.append([point_list[point_size-1], point_list[0]])

        #print(edge)

        return edge

    def calc_centroid(self):
        A = [self.vertices[i] for i in sorted(self.pints_id)]
        c = [0.0 for i in range(self.dim)]

        for i in range(self.dim):
            for j in range(self.dim):
                c[j] += A[i][j]

        for i in range(self.dim):
            c[i] /= self.dim

        # print("centroidal: {}".format(c))

        return c

vertices = [[1, 1], [3, -1]]
# vertices = [[1, 1, 2], [0, -2, 1], [3, -1, 0]]
# vertices = [[1, 1, 2, 3], [0, -2, 1, 5], [3, -1, 0, 9], [0, 0, 1, 2], [2, 3, 4, 2]]

face1 = Facet([0, 1], vertices)
# face1 = Facet([0, 1, 2], vertices)
# face1 = Facet([1, 2, 3, 4], vertices)
# face2 = Facet([0, 1, 3], vertices)
# face3 = Facet([0, 2, 3], vertices)
# face4 = Facet([1, 2, 3], vertices)

print(face1.normal)


faces = [face1]






fig = plt.figure(figsize=plt.figaspect(1.0))
ax = Axes3D(fig)


### 頂点表示
for v in vertices:
    if face1.dim == 3:
        ax.scatter(v[0], v[1], v[2], marker="o", c='r', s=60)
    if face1.dim == 2:
        ax.scatter(v[0], v[1], marker="o", c='r', s=60)

edge = []
for f in faces:
    for e in f.edge:
        edge.append(vertices[e[0]])
        edge.append(vertices[e[1]])

if face1.dim == 3:
    line_x = [e[0] for e in edge]
    line_y = [e[1] for e in edge]
    line_z = [e[2] for e in edge]
    ax.plot(line_x, line_y, line_z, "-", color='k')

if face1.dim == 2:
    line_x = [e[0] for e in edge]
    line_y = [e[1] for e in edge]
    ax.plot(line_x, line_y, "-", color='k')


### 法線表示
offset = face1.centroid

if face1.dim == 3:
    ax.quiver(offset[0], offset[1], offset[2], face1.normal[0], face1.normal[1], face1.normal[2])
if face1.dim == 2:
    ax.plot([offset[0], offset[0] + face1.normal[0]], [offset[1], offset[1] + face1.normal[1]], "-", color='b')

plt.show()