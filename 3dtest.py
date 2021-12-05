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
        self.vertices = _vertices
        self.pints_id = _id_list
        self.edge = self.set_edge()
        #self.normal = self.calc_normal(_points)
        self.normal = self.calc_normal4D(_points)

        self.centroid = self.calc_centroid(_points)

    def calc_normal4D(self):
        A = np.array([self.vertices[i] for i in sorted(self.pints_id)])
        x = np.dot(np.linalg.pinv(A), np.array([1, 1, 1, 1]))
        x = x / np.linalg.norm(x)

        print(A)

        p0 = self.vertices[self.pints_id[0]]
        p1 = self.vertices[self.pints_id[1]]
        p2 = self.vertices[self.pints_id[2]]
        p3 = self.vertices[self.pints_id[3]]

        u1 = np.array(p1) - np.array(p0)
        u2 = np.array(p2) - np.array(p1)
        u3 = np.array(p3) - np.array(p2)

        print("normal(inv): {} dot: {}".format(x, np.round(np.dot(u1, x), 3)))
        print("normal(inv): {} dot: {}".format(x, np.round(np.dot(u2, x), 3)))
        print("normal(inv): {} dot: {}".format(x, np.round(np.dot(u3, x), 3)))

        return x

    # def calc_normal(self, _points):
    #     A = np.array([_points[i] for i in sorted(self.pints_id)])
    #     x = np.dot(np.linalg.pinv(A), np.array([1, 1, 1]))
    #     x = x / np.linalg.norm(x)
    #     #x = x * 4.0
    #     #x = np.round(x, 3)
        
    #     print(A)

        
    #     p_id0 = self.pints_id[0]
    #     p_id1 = self.pints_id[1]
    #     p_id2 = self.pints_id[2]
    #     u = [_points[p_id1][0]-_points[p_id0][0], _points[p_id1][1]-_points[p_id0][1], _points[p_id1][2]-_points[p_id0][2]]
    #     v = [_points[p_id2][0]-_points[p_id0][0], _points[p_id2][1]-_points[p_id0][1], _points[p_id2][2]-_points[p_id0][2]]

    #     cros = np.cross(u, v)
        

    #     x1, x2, x3, x4, d = sympy.symbols('x1 x2 x3 x4 d')
    #     mass = sympy.solve(
    #             [A[0,0] * x1 + A[0,1] * x2 + A[0,2] * x3 + d,
    #             A[1,0] * x1 + A[1,1] * x2 + A[1,2] * x3 + d,
    #             A[2,0] * x1 + A[2,1] * x2 + A[2,2] * x3 + d], [x1, x2, x3])


    #     print("normal(inv): {} dot: {}".format(x, np.round(np.dot(u, x), 3)))
    #     print("normal(crs): {} dot: {}".format(cros, np.round(np.dot(u, cros), 3)))
    #     print("equations  : {} ".format(mass))
        

    #     return x

    def set_edge(self):
        point_list = sorted(self.pints_id)
        point_size = len(point_list)
        edge = []
        for i in range(point_size-1):
            edge.append([point_list[i], point_list[i+1]])
        edge.append([point_list[point_size-1], point_list[0]])

        return edge

    def calc_centroid(self, _points):
        A = [_points[i] for i in sorted(self.pints_id)]
        x = (A[0][0] +A[1][0]+A[2][0]) / 3.0
        y = (A[0][1] +A[1][1]+A[2][1]) / 3.0
        z = (A[0][2] +A[1][2]+A[2][2]) / 3.0

        return [x, y, z]

vertices = [[1, 1, 2, 3], [0, -2, 1, 5], [3, -1, 0, 9], [0, 0, 1, 2], [2, 3, 4, 2]]

face1 = Facet([1, 2, 3, 4], vertices)
# face2 = Facet([0, 1, 3], vertices)
# face3 = Facet([0, 2, 3], vertices)
# face4 = Facet([1, 2, 3], vertices)

print(face1.normal)


faces = [face1]






fig = plt.figure(figsize=plt.figaspect(1.0))
ax = Axes3D(fig)



for v in vertices:
    ax.scatter(v[0], v[1], v[2], marker="o", c='r', s=60)

edge = []
for f in faces:
    for e in f.edge:
        edge.append(vertices[e[0]])
        edge.append(vertices[e[1]])

line_x = [e[0] for e in edge]
line_y = [e[1] for e in edge]
line_z = [e[2] for e in edge]

ax.plot(line_x, line_y, line_z, "-", color='k')



# offset = vertices[face1.pints_id[0]]
offset = face1.centroid


ax.scatter(face1.normal[0]+offset[0], face1.normal[1]+offset[1], face1.normal[2]+offset[2], marker="o", c='b', s=60)

# ax.plot([offset[0], 2+offset[0]], [offset[1], 8+offset[1]], [offset[2], 1+offset[2]], "-", color='b') # sympyにより連立方程式の解の結果

# ax.plot([offset[0], face1.normal[0]+offset[0]], [offset[1], face1.normal[1]+offset[1]], [offset[2], face1.normal[2]+offset[2]], "-", color='y') # 逆行列による連立方程式の解の結果

ax.quiver(offset[0], offset[1], offset[2], face1.normal[0], face1.normal[1], face1.normal[2])

plt.show()