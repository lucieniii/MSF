import os
os.chdir(os.path.split(__file__)[0])

import utils
import laplacian
import optimize
import solver
import voronoi

import trimesh
import numpy as np
from tqdm import tqdm

class MSF:

    config = None
    logger = None
    mesh_saver = None
    configger = None
    laplacian = None
    solver = None
    
    mesh = None
    is_fixed = None
    voronoi_poles = None
    vertices_color = None

    scale = None
    scale_min = None
    scale_delta = None
    scale_max = None
    scale_fix = None

    alpha = None

    WL = None
    WH = None
    WM = None
    
    it_counter = 0

    __initialized = None

    def __init__(self):
        self.configger = utils.Config()
        self.config = self.configger.load_config()
        self.logger = utils.Logger()
        self.mesh_saver = utils.MeshSaver(self.logger.result_path)
        self.laplacian = laplacian.Laplacian(self.config['Laplacian_type'])
        self.solver = solver.Solver(self.config['lsqr_args'])
        self.configger.save_config(self.config, self.logger.result_path)
        self.__initialized = False

    def load_mesh(self, path):

        self.logger.logging('Loading mesh.', to_console=True)

        self.mesh = utils.read_mesh_trimesh(path)

        # Check if the mesh is watertight
        assert self.mesh.is_watertight, 'mesh is not watertight'
        
        # init is_fixed and vertices_color
        self.is_fixed = np.zeros((self.mesh.vertices.shape[0]), dtype=bool)
        self.vertices_color = np.ones_like(self.mesh.vertices)

        # Calculate scale
        bbox = self.mesh.bounding_box
        diag_len = np.max(np.linalg.norm(bbox.vertices - bbox.vertices[0], axis=1))
        self.scale = diag_len * self.config['scale']
        self.scale_min = diag_len * self.config['scale_min']
        self.scale_delta = diag_len * self.config['scale_delta']
        self.scale_max = diag_len * self.config['scale_max']
        self.scale_fix = diag_len * self.config['scale_fix']

        # Calculate Voronoi Points
        if self.config['use_reconstruction_in_voronoi']:
            self.logger.logging('Performing mesh reconstruction.', to_console=True)
            self.voronoi_poles = voronoi.calculate_voronoi_poles(
                self.laplacian.mesh_reconstruction(
                    self.mesh, self.config['k']))
        else:
            self.voronoi_poles = voronoi.calculate_voronoi_poles(self.mesh)
        voronoi_mesh = trimesh.Trimesh(vertices=self.voronoi_poles, faces=self.mesh.faces)
        self.mesh_saver.save_voronoi(voronoi_mesh)

        # set parameters
        self.WL = np.ones((self.mesh.vertices.shape[0])) * self.config['wL']
        self.WH = np.ones((self.mesh.vertices.shape[0])) * self.config['wH']
        self.WM = np.ones((self.mesh.vertices.shape[0])) * self.config['wM']

        self.alpha = self.config['alpha'] * np.pi / 180.

        # optimize the mesh one time at first
        self.optimize_mesh()

        self.__initialized = True

    def optimize_mesh(self):
    
        self.voronoi_poles, self.WL, self.WH, self.WM, self.is_fixed, self.vertices_color, n_collapsed, n_it = optimize.collapse_Edges(
            self.mesh, self.scale, self.voronoi_poles, self.WL, self.WH, self.WM, self.is_fixed, self.vertices_color)
        self.logger.logging(n_collapsed, "vertices collapsed in", n_it, "iterations.")
        
        n_added, self.voronoi_poles, n_it = optimize.split_Faces(
            self.mesh, self.alpha, self.voronoi_poles, self.is_fixed, self.scale_fix)
        self.logger.logging(n_added, "vertices added in", n_it, "iterations.")

        self.WL = np.hstack((self.WL, np.ones((n_added,)) * self.config['wL']))
        self.WH = np.hstack((self.WH, np.ones((n_added,)) * self.config['wH']))
        self.WM = np.hstack((self.WM, np.zeros((n_added,))))
        self.is_fixed = np.hstack((self.is_fixed, np.zeros((n_added,), dtype=bool)))
        for _ in range(n_added):
            self.vertices_color = np.vstack((self.vertices_color, np.array([0., 0., 1.])))

    def detect_skeleton_vertices(self):

        vertices = np.asarray(self.mesh.vertices)
        neighbours = self.mesh.vertex_neighbors
        fix_counter = 0

        for v in range(vertices.shape[0]):
            if self.is_fixed[v]:
                continue
            badCounter = 0
            for nv in neighbours[v]:
                el = np.linalg.norm(vertices[v] - vertices[nv])
                if el < self.scale_fix and not optimize.can_be_collapsed(self.mesh.vertex_neighbors, v, nv):
                    badCounter += 1
            if badCounter >= 2:
                self.is_fixed[v] = True
                fix_counter += 1

        self.WL[self.is_fixed] = 0
        self.WH[self.is_fixed] = 0
        self.WM[self.is_fixed] = 0
        self.vertices_color[self.is_fixed] = np.array([1., 0., 0.])

        self.logger.logging(fix_counter, 'vertices fixed in this interation.', to_console=True)

    def iterate(self):
        
        self.logger.logging('-------------------------------------', to_console=True)
        
        assert self.__initialized, 'Mesh not be loaded!'

        self.it_counter += 1

        self.logger.logging('Iteration %d:' % self.it_counter, to_console=True)

        self.logger.logging('Constructing Laplace matrix.')
        laplace_csr = self.laplacian.generator(self.mesh, self.is_fixed)
        
        self.logger.logging('Solving linear system.', to_console=True)
        vertices = np.asarray(self.mesh.vertices.copy())
        self.mesh.vertices = self.solver.update_parameter(
            self.WL, self.WH, self.WM, vertices, laplace_csr, self.voronoi_poles)
        
        self.mesh.vertices = self.solver.solve(vertices, logfunc=self.logger.logging)

        self.logger.logging('Perform Implicit Laplacian mesh smoothing.', to_console=True)
        self.mesh = self.laplacian.implicit_Laplacian_mesh_smoothing(
            self.mesh, self.config['smooth_lam'], self.config['smooth_it'], self.is_fixed)
        
        self.logger.logging('Optimizing mesh.', to_console=True)
        self.optimize_mesh()

        self.logger.logging('Trying to fix vertices.', to_console=True)
        self.detect_skeleton_vertices()

        if self.config['use_dynamic_scale']:
            self.scale += self.scale_delta
            if self.scale > self.scale_max:
                self.scale = self.scale_min
            self.logger.logging('Adjust scale to', self.scale, to_console=True)

        n_vertices = self.mesh.vertices.shape[0]
        n_fixed = np.where(self.is_fixed)[0].shape[0]
        self.logger.logging(n_fixed, 'of', n_vertices, 'fixed in this iteration.', to_console=True)

        voronoi_mesh = trimesh.Trimesh(vertices=self.voronoi_poles, faces=self.mesh.faces)
        self.mesh_saver.save_voronoi(voronoi_mesh)
        self.mesh_saver.save_skeleton(self.mesh, self.vertices_color)

if __name__ == '__main__':

    msf = MSF()
    msf.load_mesh("../models/armadillo.obj")
    