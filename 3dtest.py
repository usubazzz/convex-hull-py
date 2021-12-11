import re
import numpy as np
import numpy.random as rd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sympy

import copy

import facet


###     # 法線ベクトルは原点を中心としてる
### 法線ベクトルの向きは，重心よりも外を向くようにする
### そのため，simplexの重心を知る必要がある





class QuickHull():
    def __init__(self, _points, _dim):
        self.dim = _dim
        self.points = _points
        self.facets = []

        self.centroid = [0.0 for i in range(self.dim)]

        # self.asiign_facets_to_points()


    def calc_centroid(self):
        _simplex_points_id = []
        _centroid = np.array([0.0 for i in range(self.dim)])
        for f in self.facets:
            _simplex_points_id.extend(f.points_id)

        _uniq_id = list(set(_simplex_points_id))

        for _id in _uniq_id:
            _centroid += np.array(self.points[_id])
        _centroid /= len(_uniq_id)

        return _centroid

    def generate_first_simplex(self):
        _facet = self.generate_first_facet()
        print(_facet.normal)
        self.asiign_facets_to_points([_facet])
        print(_facet.out_points_id)
        farthest_point_id = self.pick_farthest_point_id(_facet)

        _facets = self.create_facets(_facet.ridge, farthest_point_id)

        self.facets.append(_facet)
        self.facets.extend(_facets)

        print(len(self.facets))

    def generate_first_facet(self):
        _maxs = [0 for i in range(self.dim)]
        _ids  = [0 for i in range(self.dim)]

        for i in range(len(self.points)):
            for j in range(len(self.points[0])):
                p_abs = np.abs(self.points[i][j])
                if p_abs > _maxs[j]:
                    _maxs[j] = p_abs
                    _ids[j] = i

        # 各次元の最大値となる点をシンプレックスの点とする
        return(facet.Facet(_ids, self.points))

    # ファセットよりも上になる点を割り当てる
    def asiign_facets_to_points(self, facets):
        for i, p in enumerate(self.points):
            for j, f in enumerate(facets):
                if i in f.points_id:
                    continue                    # ファセットを構成する点は無視
                v = np.array(p) - np.array(self.centroid)
                dot = np.dot(f.normal, v)
                if dot >= 0:
                    print("p{} -> f{}".format(i, j))
                    f.out_points_id.append(i)   # p を fに割り当てる
                    break

    # 最遠の点を取り出す
    def pick_farthest_point_id(self, facet):
        if facet.out_points_id == []:
            print("set out points")
            return None

        _max_d = 0
        _max_id = 0

        dO = np.dot(facet.normal, self.points[facet.points_id[0]])    # 原点から超平面までの距離
        for id in facet.out_points_id:
            dOP = np.dot(facet.normal, self.points[id])           # 超平面の法線に沿った原点から点Pまでの距離
            dP = np.abs(dOP - dO)                               # 超平面から点Pまでの距離
            if dP > _max_d:
                _max_d = dP
                _max_id = id

        print("d_max id: {}".format(_max_id))

        return _max_id

    # 最遠の点が超平面の上になるファセットを取り出す
    def pick_visible_facets(self, _point_id):
        _point = self.points[_point_id]
        _visible_facets = []
        for i, f in enumerate(self.facets):
            v = np.array(_point) - np.array(self.centroid)
            dot = np.dot(f.normal, v)
            # print("v: {}, n: {}, nv: {}".format(point, f.normal, dot))
            if dot > 0:
                print("visible: facet {}".format(i))
                _visible_facets.append(f)

        return _visible_facets

    # 各ファセットで共通しない点を取得(隣接するファセットの端の点)
    def get_unique_points(self, facets):
        _all_points_id = []

        for f in facets:
            _all_points_id.extend(f.points_id)

        return [x for x in set(_all_points_id) if _all_points_id.count(x) == 1]

    def create_facets(self, _ridge, _point_id):
        fs = []

        for r in _ridge:
            _points_id = []
            _points_id.append(r)
            _points_id.append(_point_id)

            f = facet.Facet(_points_id, self.points)
            fs.append(f)

        return fs


    def run(self):
        self.generate_first_simplex()


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

facets_stack_id = []
visible_set = []

points = [[(rnd_max-rnd_min)*rd.random() + rnd_min for i in range(D)] for j in range(N)]


body = QuickHull(points, 2)

body.run()

print(len(body.facets))

# facets = mesh.facets

# ## 外側の点があるファセットをスタックさせる
# if mesh.facets[0] != []:
#     facets_stack = mesh.facets

# while facets_stack !=[]:
#     f = facets_stack.pop()
#     max_id = mesh.pick_farthest_point(f, points)

#     visible_facets = mesh.pick_visible_facet(points[max_id])

#     if visible_facets != []:
#         uni = mesh.get_unique_points(visible_facets)
#         new_edge = [[max_id, uni[i]] for i in range(len(uni))]

#         for edge in new_edge:
#             new_face = facet.Facet(edge, points)

#             # new_facet に外側の点を割り当てる



#             mesh.facets.append(new_face)







fig = plt.figure(figsize=plt.figaspect(1.0))
if body.dim == 3:
    ax = Axes3D(fig)
if body.dim == 2:
    ax = fig.add_subplot(1,1,1) # 2D


### 頂点表示
for v in points:
    if body.dim == 3:
        ax.scatter(v[0], v[1], v[2], marker="o", c='r', s=60)
    if body.dim == 2:
        ax.scatter(v[0], v[1], marker="o", c='r', s=20)

# edge = []
# for f in facets:
#     for e in f.edge:
#         edge.append(points[e[0]])
#         edge.append(points[e[1]])

# if body.dim == 3:
#     line_x = [e[0] for e in edge]
#     line_y = [e[1] for e in edge]
#     line_z = [e[2] for e in edge]
#     ax.plot(line_x, line_y, line_z, "-", color='k')

if body.dim == 2:
    for f in body.facets:
        line_x = [body.points[pid][0] for pid in f.ridge]
        line_y = [body.points[pid][1] for pid in f.ridge]
        ax.plot(line_x, line_y, "-", color='k')


# ### 法線表示
# offset = body.facets[0].centroid
# face = body.facets[0]

# if body.dim == 3:
#     ax.quiver(offset[0], offset[1], offset[2], face.normal[0], face.normal[1], face.normal[2])
# if body.dim == 2:
#     for f in facets:
#         offset = f.centroid
#         ax.plot([offset[0], offset[0] + f.normal[0]], [offset[1], offset[1] + f.normal[1]], "-", color='b')


# ### 最も遠い頂点
# if body.dim == 2:
#     ax.scatter(points[max_id][0], points[max_id][1], marker="o", c='b', s=60)


plt.show()