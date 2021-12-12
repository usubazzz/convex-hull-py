import numpy as np

class Facet():
    # 法線方向の正負のためにconvexからcentroid情報がほしい?

    def __init__(self, _id_list, _points):
        self.dim = len(_id_list)
        self.points_id = _id_list

        self.ridge = []
        self.normal = self.calc_normal(_points)

        self.out_points_id = []

        self.set_ridge()

    # ToDo: 2D 以外にも対応させたい
    def calc_normal(self, _points):
        A = np.array([_points[i] for i in sorted(self.points_id)])
        x = np.dot(np.linalg.inv(A), np.ones(self.dim))
        x = x / np.linalg.norm(x)
 
        p0 = _points[self.points_id[0]]
        p1 = _points[self.points_id[1]]
        # p2 = self.points[self.points_id[2]]
        # p3 = self.points[self.pints_id[3]]

        u1 = np.array(p1) - np.array(p0)

        # print("normal(inv): {} dot: {}".format(x, np.round(np.dot(u1, x), 3)))

        return x

    def set_ridge(self):
        _dim = self.dim
        _ridge = []

        for i in range(_dim):
            for j in range(_dim-1):
                _ridge.append(self.points_id[(i+j)%_dim])
        
        self.ridge = list(_ridge)
