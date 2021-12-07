import re
import numpy as np
import numpy.random as rd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sympy


###     # 法線ベクトルは原点を中心としてる
### 法線ベクトルの向きは，重心よりも外を向くようにする
### そのため，simplexの重心を知る必要がある

class Point():
    def __init__(self):
        self.tag = None # 割り当てられてるファセット

class Facet():
    # ファセットは，頂点のidを知る必要はないかも
    def __init__(self, _id_list, _points):
        self.dim = len(_id_list)
        self.points = _points
        self.points_id = _id_list
        self.edge = self.set_edge()
        self.normal = self.calc_normal()
        self.centroid = self.calc_centroid()

        self.out_points_id = []

    def calc_normal(self):
        A = np.array([self.points[i] for i in sorted(self.points_id)])
        x = np.dot(np.linalg.inv(A), np.ones(self.dim))         # pinv である必要はない？
        x = x / np.linalg.norm(x)
 
        p0 = self.points[self.points_id[0]]
        p1 = self.points[self.points_id[1]]
        # p2 = self.points[self.points_id[2]]
        # p3 = self.points[self.pints_id[3]]

        u1 = np.array(p1) - np.array(p0)
        # u2 = np.array(p2) - np.array(p1)
        # u3 = np.array(p3) - np.array(p2)

        # ベクトルから求める方法
        #x = [-u1[1], u1[0]]

        print("normal(inv): {} dot: {}".format(x, np.round(np.dot(u1, x), 3)))
        # print("normal(inv): {} dot: {}".format(x, np.round(np.dot(u2, x), 3)))
        # print("normal(inv): {} dot: {}".format(x, np.round(np.dot(u3, x), 3)))

        return x

    # visibleから作るときのため，指定できるように改良する
    def set_edge(self):
        point_list = sorted(self.points_id)
        point_size = len(point_list)
        edge = []
        for i in range(point_size-1):
            edge.append([point_list[i], point_list[i+1]])
        edge.append([point_list[point_size-1], point_list[0]])

        #print(edge)

        return edge

    def calc_centroid(self):
        A = [self.points[i] for i in sorted(self.points_id)]
        c = [0.0 for i in range(self.dim)]

        for i in range(self.dim):
            for j in range(self.dim):
                c[j] += A[i][j]

        for i in range(self.dim):
            c[i] /= self.dim

        # print("centroidal: {}".format(c))

        return c


def pick_farthest_point(facet, _points):
    if facet.out_points_id == []:
        print("set out points")
        return None

    _max_d = 0
    _max_id = 0

    dO = np.dot(facet.normal, _points[facet.points_id[0]])    # 原点から超平面までの距離
    for id in facet.out_points_id:
        dOP = np.dot(facet.normal, _points[id])           # 超平面の法線に沿った原点から点Pまでの距離
        dP = np.abs(dOP - dO)                               # 超平面から点Pまでの距離
        if dP > _max_d:
            _max_d = dP
            _max_id = id

    print("d_max id: {}".format(_max_id))

    return _max_id

# 最遠の点が超平面の上になるファセットを取り出す
def pick_visible_facet(facets, point):

    
    _visible_facets = []
    for i, f in enumerate(facets):
        v = np.array(point) - np.array(f.centroid)
        dot = np.dot(f.normal, v)
        # print("v: {}, n: {}, nv: {}".format(point, f.normal, dot))
        if dot > 0:
            print("visible: facet {}".format(i))
            _visible_facets.append(f)

    return _visible_facets

# 各ファセットで共通しない点を取得(隣接するファセットの端の点)
def get_unique_points(facets):
    _all_points_id = []

    for f in facets:
        _all_points_id.extend(f.points_id)

    return [x for x in set(_all_points_id) if _all_points_id.count(x) == 1]



# visibleファセットを削除，ユニークな点と最遠の点から新たなファセットを作成してプッシュする
def create_facets(_points, edge_points_id):

    f = Facet(edge_points_id, _points)



N = 10
D = 2

# ランダム範囲 min <= rnd < max
rnd_min = -10
rnd_max = 10
seed = rd.randint(0, 300)
# rd.seed(seed)
# rd.seed(230)
rd.seed(102)
print("seed: {}".format(seed))

points = [[(rnd_max-rnd_min)*rd.random() + rnd_min for i in range(D)] for j in range(N)]

# points = [[1, 1], [3, -1], [2, 4], [4, 5], [-3, -5]]
# points = [[1, 1, 2], [0, -2, 1], [3, -1, 0]]
# points = [[1, 1, 2, 3], [0, -2, 1, 5], [3, -1, 0, 9], [0, 0, 1, 2], [2, 3, 4, 2]]

face1 = Facet([0, 1], points)
face2 = Facet([1, 2], points)
# face1 = Facet([0, 1, 2], points)
# face1 = Facet([1, 2, 3, 4], points)
# face2 = Facet([0, 1, 3], points)
# face3 = Facet([0, 2, 3], points)
# face4 = Facet([1, 2, 3], points)


facets = [face1, face2]


# 各頂点の割当
for i, p in enumerate(points):
    for j, f in enumerate(facets):
        if i in f.points_id:
            # ファセットを構成する点は無視
            continue
        v = np.array(p) - np.array(f.centroid)
        dot = np.dot(f.normal, v)
        if dot >= 0:
            # p を fに割り当てる
            print("p{} -> f{}".format(i, j))
            f.out_points_id.append(i)
            break

max_id = pick_farthest_point(face1, points)


visible_facets = pick_visible_facet(facets, points[max_id])


g = get_unique_points(visible_facets)
print(g)


in_points = []    # convex の内側にある頂点



fig = plt.figure(figsize=plt.figaspect(1.0))
if face1.dim == 3:
    ax = Axes3D(fig)
if face1.dim == 2:
    ax = fig.add_subplot(1,1,1) # 2D


### 頂点表示
for v in points:
    if face1.dim == 3:
        ax.scatter(v[0], v[1], v[2], marker="o", c='r', s=60)
    if face1.dim == 2:
        ax.scatter(v[0], v[1], marker="o", c='r', s=20)

edge = []
for f in facets:
    for e in f.edge:
        edge.append(points[e[0]])
        edge.append(points[e[1]])

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
    for f in facets:
        offset = f.centroid
        ax.plot([offset[0], offset[0] + f.normal[0]], [offset[1], offset[1] + f.normal[1]], "-", color='b')


### 最も遠い頂点
if face1.dim == 2:
    ax.scatter(points[max_id][0], points[max_id][1], marker="o", c='b', s=60)


plt.show()