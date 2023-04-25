from scipy.sparse.linalg import lsqr
from scipy.sparse import vstack

import numpy as np

class Solver:

    WL = None
    WH = None
    WM = None
    L = None
    vertices = None
    voronoi_poles = None
    lsqr_args = None

    def __init__(self, lsqr_args):
        self.lsqr_args = lsqr_args

    def update_parameter(self, WL, WH, WM, vertices, L, voronoi_poles):
        self.WL = WL
        self.WH = WH
        self.WM = WM
        self.L = L
        self.vertices = vertices
        self.voronoi_poles = voronoi_poles

    def __construct_equation(self):
        vn = self.vertices.shape[0]
        A = vstack((self.WL * self.L, self.WH, self.WM))
        b = np.vstack((np.zeros((vn, 3)), self.WH * self.vertices, self.WM * self.voronoi_poles))
        return A, b

    def solve(self, x0, logfunc=None, log_to_console=False):
        
        A, b = self.__construct_equation()

        new_vertices0 = lsqr(A=A, b=b[:, 0], **self.lsqr_args, x0=x0[:, 0])
        new_vertices1 = lsqr(A=A, b=b[:, 1], **self.lsqr_args, x0=x0[:, 1])
        new_vertices2 = lsqr(A=A, b=b[:, 2], **self.lsqr_args, x0=x0[:, 2])

        new_vertices = np.zeros(self.vertices.shape)
        new_vertices[:, 0] = new_vertices0[0]
        new_vertices[:, 1] = new_vertices1[0]
        new_vertices[:, 2] = new_vertices2[0]

        if logfunc:
            logfunc('lsqr status:', new_vertices0[1], 'l2 norm', np.linalg.norm(A @ new_vertices0[0] - b[:, 0]), to_console=log_to_console)
            logfunc('lsqr status:', new_vertices1[1], 'l2 norm', np.linalg.norm(A @ new_vertices1[0] - b[:, 1]), to_console=log_to_console)
            logfunc('lsqr status:', new_vertices2[1], 'l2 norm', np.linalg.norm(A @ new_vertices2[0] - b[:, 2]), to_console=log_to_console)

        return new_vertices