from scipy.spatial import Voronoi
import numpy as np
from tqdm import tqdm

INFINITY = 1e9

def calculate_voronoi_poles(mesh):
    
    vertices = np.asarray(mesh.vertices.copy())
    vertice_normals = np.asarray(mesh.vertex_normals.copy())

    print("Calculating Voronoi Diagram.")
    vor = Voronoi(vertices)
    vor_centers = vor.vertices
    cells = vor.regions
    cell_indices = vor.point_region

    vertices_total = vertices.shape[0]
    voronoi_poles = np.zeros(vertices.shape)

    for vi in tqdm(range(vertices_total), desc='Set Voronoi pole for each vertex'):
        vor_cell = cells[cell_indices[vi]]
        vertice = vertices[vi]
        vertice_normal = vertice_normals[vi]
        max_neg_proj = INFINITY
        voronoi_pole = None
        for vci in vor_cell:
            vor_center = vor_centers[vci]
            if vci == -1:
                continue
            proj = np.dot(vor_center - vertice, vertice_normal)
            if proj < max_neg_proj:
                max_neg_proj = proj
                voronoi_pole = vor_center
        voronoi_poles[vi] = voronoi_pole

    return voronoi_poles