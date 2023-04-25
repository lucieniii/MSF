import numpy as np
import trimesh
from scipy.sparse.linalg import eigs, spsolve
from scipy.sparse import csr_matrix, coo_matrix, identity
from tqdm import tqdm
from enum import Enum

class Laplacian:

    EPS = 1e-8
    generator = None

    def __init__(self, laplacian_type):
        super().__init__()
        if laplacian_type == 'uniform':
            self.generator = self.trimesh_generate_uniform_Laplace_matrix
        elif laplacian_type == 'tangent':
            self.generator = self.trimesh_generate_tangent_Laplace_matrix
        elif laplacian_type == 'cotangent':
            self.generator = self.trimesh_generate_cotangent_Laplace_matrix

    def trimesh_generate_area_list(self, mesh):

        area_list = np.zeros((len(mesh.vertices)))
        face_count = mesh.faces.shape[0]
        face_area = np.asarray(mesh.area_faces)
        faces = np.asarray(mesh.faces)
        for i in range(face_count):
            area_list[faces[i]] += face_area[i] / 3
        
        return area_list

    def trimesh_generate_cotangent_Laplace_matrix(self, mesh, is_fixed):
        
        vertices_face_indexs = [[0,1,2],[1,0,2],[2,0,1]] 
        laplace_dict = {}
        area_list = self.trimesh_generate_area_list(mesh)
        face_angles = np.asarray(mesh.face_angles)
        faces = np.asarray(mesh.faces)
        face_count = faces.shape[0]
        
        with tqdm(total=face_count) as tbar:
            tbar.set_description('Constructing Cotangent Laplace Matrix.')
            for face, angles in zip(faces, face_angles):
                for i in range(3):
                    current_angle = angles[i]
                    v0_index, v1_index, v2_index = face[vertices_face_indexs[i]]
                    delta = 1 / (area_list[v0_index] * np.tan(current_angle) + self.EPS)
                    if not is_fixed[v1_index]:
                        laplace_dict[(v1_index, v1_index)] = laplace_dict.get((v1_index, v1_index), 0) - delta
                        laplace_dict[(v1_index, v2_index)] = laplace_dict.get((v1_index, v2_index), 0) + delta
                    if not is_fixed[v2_index]:
                        laplace_dict[(v2_index, v2_index)] = laplace_dict.get((v2_index, v2_index), 0) - delta
                        laplace_dict[(v2_index, v1_index)] = laplace_dict.get((v2_index, v1_index), 0) + delta
                tbar.update(1)

        # Construct CSR Matrix
        rows, cols = zip(*laplace_dict.keys())
        values = list(laplace_dict.values())
        coo = coo_matrix((values, (rows, cols)), shape=(len(mesh.vertices), len(mesh.vertices)))
        csr = coo.tocsr()

        return csr

    def trimesh_generate_tangent_Laplace_matrix(self, mesh, is_fixed):
        
        vertices_face_indexs = [[0,1,2],[1,0,2],[2,0,1]] 
        laplace_dict = {}
        face_angles = np.asarray(mesh.face_angles)
        faces = np.asarray(mesh.faces)
        face_count = faces.shape[0]
        vertices = np.asarray(mesh.vertices)
        
        with tqdm(total=face_count) as tbar:
            tbar.set_description('Constructing Tangent Laplace Matrix.')
            for face, angles in zip(faces, face_angles):
                for i in range(3):
                    alpha0, alpha1, alpha2 = angles[vertices_face_indexs[i]]
                    v0_index, v1_index, v2_index = face[vertices_face_indexs[i]]
                    el = np.linalg.norm(vertices[v1_index] - vertices[v2_index]) + self.EPS
                    if not is_fixed[v1_index]:
                        laplace_dict[(v1_index, v1_index)] = laplace_dict.get((v1_index, v1_index), 0) - np.tan(alpha1 / 2.) / el
                        laplace_dict[(v1_index, v2_index)] = laplace_dict.get((v1_index, v2_index), 0) + np.tan(alpha1 / 2.) / el
                    if not is_fixed[v2_index]:
                        laplace_dict[(v2_index, v2_index)] = laplace_dict.get((v2_index, v2_index), 0) - np.tan(alpha2 / 2.) / el
                        laplace_dict[(v2_index, v1_index)] = laplace_dict.get((v2_index, v1_index), 0) + np.tan(alpha2 / 2.) / el
                tbar.update(1)

        # Construct CSR Matrix
        rows, cols = zip(*laplace_dict.keys())
        values = list(laplace_dict.values())
        coo = coo_matrix((values, (rows, cols)), shape=(len(mesh.vertices), len(mesh.vertices)))
        csr = coo.tocsr()

        return csr
    
    def trimesh_generate_uniform_Laplace_matrix(self, mesh, is_fixed):
        """
        Compute the matrix of uniform Laplace-Beltrami operator for a mesh.

        Args:
            mesh:       Trimesh, the data of mesh loaded by Trimesh.
        
        Returns:
            csr_L:      scipy.sparse.csr_matrix, 
                        the sparse representation of the uniform Laplace-Beltrami matrix 
        """
        
        laplace_dict = {}
        neighbours = mesh.vertex_neighbors
        with tqdm(total=len(neighbours)) as tbar:
            tbar.set_description('Constructing Uniform Laplace Matrix.')
            for i, n in enumerate(neighbours):
                if is_fixed[i]:
                    continue
                nn = len(n)
                laplace_dict[(i, i)] = -1.
                for ni in n:
                    laplace_dict[(i, ni)] = 1. / nn

                tbar.update(1)
        
        # Construct CSR Matrix
        rows, cols = zip(*laplace_dict.keys())
        values = list(laplace_dict.values())
        coo = coo_matrix((values, (rows, cols)), shape=(len(mesh.vertices), len(mesh.vertices)))
        csr = coo.tocsr()

        return csr
    
    def implicit_Laplacian_mesh_smoothing(self, mesh, lam=0.01, t=20, is_fixed=None):
        """
        Implicit Laplacian mesh smoothing.

        Args:
            mesh:           Trimesh, the data of mesh loaded by Trimesh.
            lam:            float, represents lambda
            t:              int, iteration round
        
        Returns:
            smooth_mesh:    Trimesh, smoothed mesh.
        """
        
        # Copy the original mesh
        smooth_mesh = mesh.copy()
        
        # Get the vertices corrdinates
        smooth_vertices = np.asarray(mesh.vertices)
        vn = smooth_vertices.shape[0]

        # Compute Laplace-Beltrami matrix
        csr_uL = self.generator(mesh, is_fixed)

        # Iterate to smooth the mesh.
        for _ in tqdm(range(t), desc='Implicit Laplacian mesh smoothing'):
            # Solve the linear system
            smooth_vertices = spsolve(identity(vn) - lam * csr_uL, smooth_vertices)
            smooth_mesh.vertices = smooth_vertices

        return smooth_mesh
    
    def __construct_Laplace_Beltrami_matrix(self, mesh):
        """
        Compute the matrix of un-uniform Laplace-Beltrami operator for a mesh.

        Args:
            mesh:       Trimesh, the data of mesh loaded by Trimesh.
        
        Returns:
            csr_L:      scipy.sparse.csr_matrix, 
                        the sparse representation of the un-uniform Laplace-Beltrami matrix 
        """

        # To contruct csr_matrix
        data = []
        indices = []
        indptr = [0]

        vertices = mesh.vertices
        neighbours = mesh.vertex_neighbors
        for i, n in enumerate(neighbours):
            n_cnt = len(n)
            rdata = np.ones((n_cnt + 1,))
            rdata[n_cnt] = -n_cnt
            rindices = np.zeros((n_cnt + 1,))
            rindices[:n_cnt] = np.array(n)
            rindices[n_cnt] = i
            ri = np.argsort(rindices)
            data.extend(rdata[ri].tolist())
            indices.extend(rindices[ri].tolist())
            indptr.append(indptr[-1] + n_cnt + 1)
            
        csr_L = csr_matrix((data, indices, indptr), shape=(vertices.shape[0], vertices.shape[0]), dtype=np.float32)

        return csr_L
    
    def mesh_reconstruction(self, mesh, k):

        # Compute sparse Laplace-Beltrami matrix
        csr_L = self.__construct_Laplace_Beltrami_matrix(mesh)
        
        # Compute the k-smallest eigenvalues and their eigenvectors
        w, v = eigs(csr_L, k=k, which='SM')
        w = w.real
        v = v.real

        # Get the vertices corrdinates
        vertices = np.asarray(mesh.vertices)

        # To store the corrdinates after reconstruction
        smooth_vertices = np.zeros(vertices.shape)

        # Reconstrction
        ef = vertices.T @ v
        for ch in range(3):
            smooth_vertices[:, ch] = np.sum(np.tile(ef[ch], (vertices.shape[0], 1)) * v, axis=1)
        
        smooth_mesh = trimesh.Trimesh(vertices=smooth_vertices, faces=mesh.faces)
        
        return smooth_mesh