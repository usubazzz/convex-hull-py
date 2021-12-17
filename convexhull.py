#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np
import facet

### ToDo: out points を毎回すべて計算してるので，最適化したほうが速度でる

class QuickHull():
    def __init__(self, _points, _dim):
        self.dim = _dim
        self.points = _points
        self.facets = []
        self.stacks = []

        self.centroid = [0.0 for i in range(self.dim)]

        self.facet_num = 0

    def create_facet(self, _ids, _points):
        f = facet.Facet(_ids, _points, self.facet_num)
        self.facet_num += 1
        return f

    def calc_centroid(self, _facets):
        _vertexs = []
        _centroid = np.array([0.0 for i in range(self.dim)])

        #頂点リスト作成
        for f in _facets:
            for r in f.ridge:
                _vertexs.extend(r)

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
        return(self.create_facet(_ids, self.points))

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
        # print("check facets: {}".format([f.id for f in _facets]))
        for i, f in enumerate(_facets):
            v = np.array(_point) - np.array(self.points[f.points_id[0]])
            dot = np.dot(f.normal, v)
            # print("v: {}, n: {}, nv: {}".format(self.points[f.points_id[0]], f.normal, dot))
            if dot > 0:
                _visible_facets_id.append(i)

        return _visible_facets_id

    # 重複の無い(順不同)ridgeを取得
    # [0, 1, 1, 2] -> [0, 2]
    def get_only_ridges_randorder(self, facets):
        _all_ridge = []

        for f in facets:
            for r in f.ridge:
                _all_ridge.append(r)

        print("all ridge: {}".format(_all_ridge))

        double_r = []
        result = []
        for line in _all_ridge:
            if set(line) in result:
                double_r.append(set(line))
            else:
                result.append(set(line))

        _uniq_ridges = []
        for line in _all_ridge:
            if set(line) not in double_r:
                _uniq_ridges.append(line)

        print("uniq")
        print(_uniq_ridges)
        print(double_r)

        return _uniq_ridges

    # 重複しない(順不同)vertexを取得
    # [0, 1, 1, 2] -> [0, 1, 2]
    def get_unique_vertexs_randorder(self, _vertexs):
        result = []
        for x in set(_vertexs):
            result.append(x)

        return result

    # 端のridgeと最遠点から，複数のファセットを作成
    ### 2D専用になってる
    def create_facets(self, _ridge, _point_id):
        fs = []

        for r in _ridge:
            _points_id = []
            _points_id.extend(r)
            _points_id.append(_point_id)

            print("Create points: {}".format(_points_id))
            f = self.create_facet(_points_id, self.points)
            fs.append(f)

        return fs

    # 重心から外側の方向が正になるように法線を反転させる
    def calc_facet_normal_centroid(self, _facets):
        for f in _facets:
            v = np.array(self.points[f.points_id[0]]) - np.array(self.centroid)
            dot = np.dot(f.normal, v)
            if dot < 0:
                f.normal *= -1.0

    # def delete_inpoint(self, _facets):
            
    #     for i in reversed(range(len(self.points))):
    #         _ch_in = True
    #         for f in _facets:
    #             v = np.array(self.points[i]) - np.array(self.points[f.points_id[0]])
    #             dot = np.dot(f.normal, v)
    #             if dot >= 0.0:
    #                 _ch_in = False
    #         if _ch_in:
    #             print("IN point: {}".format(i))
    #             # 削除したあと，out_pointsに登録してたやつなどどうするのか
    #             # out_pointsの点のみを扱えば問題ない
    #             # del self.points[i]

    def make_hull_step(self, _facet):
        print("=======")
        print("facets id: {}".format([f.id for f in self.facets]))
        print("stacks id: {}".format([f.id for f in self.stacks]))
        
        farthest_point_id = self.pick_farthest_point_id(_facet)

        visible_check_facets = self.stacks + self.facets
        visible_ids = self.pick_visible_facets_id(farthest_point_id, visible_check_facets)
        visible_facets = [visible_check_facets[i] for i in visible_ids]

        print("Visible facets id: {}".format([visible_check_facets[i].id for i in visible_ids]))
        print("visivle ridge: {}".format([visible_check_facets[i].ridge for i in visible_ids]))

        _only_ridges = self.get_only_ridges_randorder(visible_facets)  

        _new_facets = []
        _facets = self.create_facets(_only_ridges, farthest_point_id)
        _new_facets.extend(_facets)

        self.calc_centroid(_new_facets)
        self.calc_facet_normal_centroid(_new_facets)
        self.asiign_facets_to_points(_new_facets)

        print("NEW facets id: {}".format([f.id for f in _new_facets]))

        # 新たなファセットの元を削除
        for id in sorted(visible_ids, reverse=True):
            print("DEL ID: {}".format(visible_check_facets[id].id))
            del visible_check_facets[id]


        self.facets.clear()
        self.stacks.clear()
        for f in visible_check_facets:
            if len(f.out_points_id) > 0:
                self.stacks.append(f)
            else:
                self.facets.append(f)

        for i, f in enumerate(_new_facets):
            if f.out_points_id != []:
                # print("Add stacks: {}".format(f.id))
                self.stacks.append(f)
            else:
                # print("Add facets: {}".format(f.id))
                self.facets.append(f)

        print("-------")
        print("facets id: {}".format([f.id for f in self.facets]))
        print("stacks id: {}".format([f.id for f in self.stacks]))

    def make_hull(self):
        while self.stacks != []:
            print("STACK: {}".format(len(self.stacks)))
            
            _facet = self.stacks[-1]

            self.make_hull_step(_facet)

    # 初期用のシンプレックスの作成
    def first_hull(self):
        first_facets = self.generate_first_simplex()
        self.calc_centroid(first_facets)
        self.calc_facet_normal_centroid(first_facets)
        
        self.asiign_facets_to_points(first_facets)

        
        for f in first_facets:
            if f.out_points_id != []:
                self.stacks.append(f)
            else:
                self.facets.append(f) # Convex 決定

    def run(self):
        # First
        self.first_hull()

        # Main
        self.make_hull()


