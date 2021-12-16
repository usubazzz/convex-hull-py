import numpy as np

class Facet():
    def __init__(self, _id_list, _points, id=None):
        self.dim = len(_id_list)
        self.points_id = _id_list

        self.ridge = []
        self.normal = self.calc_normal(_points)

        self.out_points_id = []

        self.id = id

        self.set_ridge()

    def calc_normal(self, _points):
        A = np.array([_points[i] for i in sorted(self.points_id)])

        # x = np.dot(np.linalg.inv(A), np.ones(self.dim))
        x = np.dot(np.linalg.pinv(A), np.ones(self.dim))
        x = x / np.linalg.norm(x)

        return x

    def set_ridge(self):
        _dim = self.dim
        _ridge = []
        
        for i in range(_dim):
            _rg = []
            for j in range(_dim-1):
                _rg.append(self.points_id[(i+j)%_dim])
            _ridge.append(_rg)
        
        self.ridge = list(_ridge)
