import re
import numpy as np
import numpy.random as rd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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
        self.stacks = []

        self.centroid = [0.0 for i in range(self.dim)]

        # self.asiign_facets_to_points()


    def calc_centroid(self, _facets):
        _vertexs = []
        _centroid = np.array([0.0 for i in range(self.dim)])

        #頂点リスト作成
        for f in _facets:
            _vertexs.extend(f.ridge)

        _uniq_vertexs = self.get_unique_vertexs_randorder(_vertexs)

        for _id in _uniq_vertexs:
            _centroid += np.array(self.points[_id])
        _centroid /= len(_uniq_vertexs)

        self.centroid = _centroid

    def generate_first_simplex(self):
        _facet = self.generate_first_facet()

        self.asiign_facets_to_points([_facet])
        print("outpoints: {}".format(_facet.out_points_id))

        # 最遠点がない場合は，法線を反転
        if _facet.out_points_id == []:
            print("flip normal")
            _facet.normal *= -1.0
            self.asiign_facets_to_points([_facet])
            print("outpoints: {}".format(_facet.out_points_id))

        farthest_point_id = self.pick_farthest_point_id(_facet)

        _facets = self.create_facets(_facet.ridge, farthest_point_id)

        _resurt_facets = []
        _resurt_facets.append(_facet)
        _resurt_facets.extend(_facets)

        return _resurt_facets

    def generate_first_facet(self):
        _maxs = [0 for i in range(self.dim)]
        _ids  = [0 for i in range(self.dim)]

        for i in range(len(self.points)):
            for j in range(len(self.points[0])):
                p_abs = np.abs(self.points[i][j])
                if p_abs > _maxs[j]:
                    if not i in _ids:
                        _maxs[j] = p_abs
                        _ids[j] = i

        # 各次元の最大値となる点をシンプレックスの点とする
        return(facet.Facet(_ids, self.points))

    # ファセットよりも上になる点を割り当てる
    def asiign_facets_to_points(self, facets):
        for f in facets:
            f.out_points_id.clear()

        for i, p in enumerate(self.points):
            for j, f in enumerate(facets):
                if i in f.points_id:
                    continue                    # ファセットを構成する点は無視
 
                v = np.array(p) - np.array(self.points[f.points_id[0]])
                
                dot = np.dot(f.normal, v)
                if dot >= 0:
                    # print("p{} -> f{}".format(i, j))
                    f.out_points_id.append(i)   # p を fに割り当てる
                    break

    # 最遠の点を取り出す
    def pick_farthest_point_id(self, facet):
        if facet.out_points_id == []:
            print("Err: No outside point")
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

        # print("d_max id: {}".format(_max_id))

        return _max_id

    # 最遠の点が超平面の上になるファセットを取り出す
    def pick_visible_facets_id(self, _point_id, _facets):
        _point = self.points[_point_id]
        _visible_facets_id = []
        for i, f in enumerate(_facets):
            v = np.array(_point) - np.array(self.points[f.points_id[0]])
            dot = np.dot(f.normal, v)
            # print("v: {}, n: {}, nv: {}".format(point, f.normal, dot))
            if dot > 0:
                print("visible: facet {}".format(i))
                _visible_facets_id.append(i)

        return _visible_facets_id

    # 重複の無い(順不同)ridgeを取得
    # [0, 1, 1, 2] -> [0, 2]
    def get_only_ridges_randorder(self, facets):
        _all_ridge = []

        for f in facets:
            _all_ridge.append(f.ridge)

        result = []
        _uniq_ridges = []
        for line in _all_ridge:
            if set(line) not in result:
                result.append(set(line))
                if _all_ridge.count(line) == 1:    # 重複しないもののみ
                    _uniq_ridges.append(line)

        return _uniq_ridges

    # 重複しない(順不同)vertexを取得
    # [0, 1, 1, 2] -> [0, 1, 2]
    def get_unique_vertexs_randorder(self, _vertexs):
        result = []
        for x in set(_vertexs):
            result.append(x)

        return result

    # 端のridgeと最遠点から，複数のファセットを作成
    def create_facets(self, _ridge, _point_id):
        fs = []

        for r in _ridge:
            _points_id = []
            _points_id.append(r)
            _points_id.append(_point_id)

            f = facet.Facet(_points_id, self.points)
            fs.append(f)

        return fs

    # 重心から外側の方向が正になるように法線を反転させる
    def calc_facet_normal_centroid(self, _facets):
        for f in _facets:
            v = np.array(self.points[f.points_id[0]]) - np.array(self.centroid)
            dot = np.dot(f.normal, v)
            if dot < 0:
                f.normal *= -1.0

    def delete_inpoint(self, _facets):
            
        for i in reversed(range(len(self.points))):
            _ch_in = True
            for f in _facets:
                v = np.array(self.points[i]) - np.array(self.points[f.points_id[0]])
                dot = np.dot(f.normal, v)
                if dot >= 0.0:
                    _ch_in = False
            if _ch_in:
                print("IN point: {}".format(i))
                # 削除したあと，out_pointsに登録してたやつなどどうするのか
                # out_pointsの点のみを扱えば問題ない
                # del self.points[i]

    def make_hull_step(self, _facet):
        farthest_point_id = self.pick_farthest_point_id(_facet)

        visible_ids = self.pick_visible_facets_id(farthest_point_id, self.stacks)
        visible_facets = [self.stacks[i] for i in visible_ids]
        _only_ridges = self.get_only_ridges_randorder(visible_facets)

        _new_facets = []
        for _ridges in _only_ridges:
            _facets = self.create_facets(_ridges, farthest_point_id)
            _new_facets.extend(_facets)

        # full_facets = []
        # full_facets.extend(self.stacks)
        # full_facets.extend(self.facets)
        # full_facets.extend(_new_facets)
        self.calc_centroid(self.stacks)
        self.calc_facet_normal_centroid(_new_facets) ## 外側の無いファセット+スタック中のファセット+新たに作ったファセット-visibleファセットに変更予定
        self.asiign_facets_to_points(_new_facets)

        for f in _new_facets:
            if f.out_points_id != []:
                self.stacks.append(f)
            else:
                self.facets.append(f)

        # 新たなファセットの元を削除
        for id in visible_ids:
            print("DEL ID: {}".format(id))
            del self.stacks[id]

    def make_hull(self):
        #### ToDo: 可視ファセットを削除するとき，スタックのどのファセットかわからない問題
        # Main
        while self.stacks != []:
            print("STACK: {}".format(len(self.stacks)))
            
            # _facet = self.stacks.pop()
            _facet = self.stacks[-1]

            self.make_hull_step(_facet)

            # farthest_point_id = self.pick_farthest_point_id(_facet)

            # visible_ids = self.pick_visible_facets_id(farthest_point_id, self.stacks)
            # visible_facets = [self.stacks[i] for i in visible_ids]
            # _only_ridges = self.get_only_ridges_randorder(visible_facets)

            # _new_facets = []
            # for _ridges in _only_ridges:
            #     _facets = self.create_facets(_ridges, farthest_point_id)
            #     _new_facets.extend(_facets)

            # # full_facets = []
            # # full_facets.extend(self.stacks)
            # # full_facets.extend(self.facets)
            # # full_facets.extend(_new_facets)
            # self.calc_centroid(self.stacks)
            # self.calc_facet_normal_centroid(_new_facets) ## 外側の無いファセット+スタック中のファセット+新たに作ったファセット-visibleファセットに変更予定
            # self.asiign_facets_to_points(_new_facets)

            # for f in _new_facets:
            #     if f.out_points_id != []:
            #         self.stacks.append(f)
            #     else:
            #         self.facets.append(f)

            # # 新たなファセットの元を削除
            # for id in visible_ids:
            #     print("DEL ID: {}".format(id))
            #     del self.stacks[id]

    def run(self):
        # First
        first_facets = self.generate_first_simplex()
        self.calc_centroid(first_facets)
        self.calc_facet_normal_centroid(first_facets)
        
        self.delete_inpoint(first_facets) # 内側の点を計算から除外する, out_pointsを使うならいらない
        
        self.asiign_facets_to_points(first_facets)

        
        for f in first_facets:
            if f.out_points_id != []:
                self.stacks.append(f)
                print("OUT P: {}".format(f.out_points_id))
            else:
                self.facets.append(f) # Convex 決定

        print("STACK NUM: {}".format(len(self.stacks)))



def plot_ridge(facets, points, ax, c):
    facets_line_x = []
    facets_line_y = []
    for f in facets:
        print("Ridge: {}".format(f.ridge))
        line_x = [points[pid][0] for pid in f.ridge]
        line_y = [points[pid][1] for pid in f.ridge]
        facets_line_x.extend(line_x)
        facets_line_y.extend(line_y)
    im = ax.plot(facets_line_x, facets_line_y, "-", color=c)
    return im

N = 30
D = 2

# ランダム範囲 min <= rnd < max
rnd_min = -10
rnd_max = 10
seed = rd.randint(0, 300)
rd.seed(seed)
# rd.seed(230)
# rd.seed(102)
# rd.seed(282)
# rd.seed(82) # Err "e:\repositories\convex-hull-py\facet.py", line 18, in calc_normal
# rd.seed(82) # Err IndexError: list assignment index out of range
print("seed: {}".format(seed))

facets_stacks_id = []
visible_set = []

points = [[(rnd_max-rnd_min)*rd.random() + rnd_min for i in range(D)] for j in range(N)]


fig = plt.figure(figsize=plt.figaspect(1.0))
ax = fig.add_subplot(1,1,1)

ims = []

body = QuickHull(points, 2)

body.run()

### 頂点表示
for v in points:
    ax.scatter(v[0], v[1], marker="o", c='r', s=20)

### 初回ファセット
im1 = plot_ridge(body.facets, body.points, ax, 'k')
im2 = plot_ridge(body.stacks, body.points, ax, 'b')
ims.append(im1+im2)



while body.stacks != []:
    _facet = body.stacks[-1]

    body.make_hull_step(_facet)

    im1 = plot_ridge(body.facets, body.points, ax, 'k')
    im2 = plot_ridge(body.stacks, body.points, ax, 'b')
    ims.append(im1+im2)


# fig = plt.figure(figsize=plt.figaspect(1.0))
# if body.dim == 3:
#     ax = Axes3D(fig)
# if body.dim == 2:
#     ax = fig.add_subplot(1,1,1) # 2D


# ### 頂点表示
# for v in points:
#     if body.dim == 3:
#         ax.scatter(v[0], v[1], v[2], marker="o", c='r', s=60)
#     if body.dim == 2:
#         ax.scatter(v[0], v[1], marker="o", c='r', s=20)


# if body.dim == 2:
#     for f in body.facets:
#         line_x = [body.points[pid][0] for pid in f.ridge]
#         line_y = [body.points[pid][1] for pid in f.ridge]
#         ax.plot(line_x, line_y, "-", color='k')


# ### 法線表示
# if body.dim == 3:
#     ax.quiver(offset[0], offset[1], offset[2], face.normal[0], face.normal[1], face.normal[2])
# if body.dim == 2:
#     for f in body.facets:
#         centroid = np.array([0.0, 0.0])
#         for pid in f.points_id:
#             centroid += np.array(body.points[pid])
#         offset = centroid / 2.0
#         ax.plot([offset[0], offset[0] + f.normal[0]], [offset[1], offset[1] + f.normal[1]], "-", color='b')

# 重心
# if body.dim == 2:
#     ax.scatter(body.centroid[0], body.centroid[1], marker="o", c='g', s=60)

# ### 最も遠い頂点
# if body.dim == 2:
#     ax.scatter(points[max_id][0], points[max_id][1], marker="o", c='b', s=60)

plt.xlim([-10, 10])
plt.ylim([-10, 10])

ani = animation.ArtistAnimation(fig, ims, interval=1000)

plt.show()