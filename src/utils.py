import os
import time
import random
import trimesh
import open3d as o3d
import numpy as np
import json

def show_mesh_o3d(plys):
    o3d.visualization.draw_geometries(plys)

def write_mesh_o3d(path, mesh):
    o3d.io.write_triangle_mesh(path, mesh)

def read_mesh_o3d(mesh_fp):
    return o3d.io.read_triangle_mesh(mesh_fp)

def read_mesh_trimesh(path, process=True):
    return trimesh.load(path, process=process)

def write_trimesh(mesh, path):
    o3d_mesh = mesh.as_open3d
    write_mesh_o3d(path, o3d_mesh)

def write_trimesh_with_color(mesh, path, colors):
    o3d_mesh = mesh.as_open3d
    o3d_mesh.vertex_colors = o3d.utility.Vector3dVector(colors)
    write_mesh_o3d(path, o3d_mesh)

class logger:

    result_path = None
    log_file = None

    def __init__(self):
        self.generate_result_path()
        self.generate_log_file()
        
    def generate_result_path(self):
        cur_dir = os.path.join(os.path.split(__file__)[0])
        result_dir = os.path.join(cur_dir, '..', 'result')
        if not os.path.exists(result_dir):
            os.mkdir(result_dir)
        dir_name = "%s__%s" % (time.strftime("%Y-%m-%d-%H_%M_%S", time.localtime()), ''.join(random.choices('ABCDEFG', k=4)))
        self.result_path = os.path.join(result_dir, dir_name)
        os.mkdir(self.result_path)

    def generate_log_file(self):

        log_file_path = os.path.join(self.result_path, 'log.txt')
        self.log_file = open(log_file_path, 'w')
    
    def logging(self, *inputs, to_console=False):
        
        print(*inputs, file=self.log_file)
        if to_console:
            print(*inputs)


class MeshSaver:

    voronoi_path = None
    skeleton_path = None
    voronoi_counter = 0
    skeleton_counter = 0
    
    def __init__(self, result_dir):
        self.voronoi_path = os.path.join(result_dir, 'voronoi')
        if not os.path.exists(self.voronoi_path):
            os.mkdir(self.voronoi_path)
        self.skeleton_path = os.path.join(result_dir, 'skeleton')
        if not os.path.exists(self.skeleton_path):
            os.mkdir(self.skeleton_path)

    def save_voronoi(self, mesh):
        file_name = 'voronoi_%d.obj' % self.voronoi_counter
        write_trimesh(mesh, os.path.join(self.voronoi_path, file_name))
        self.voronoi_counter += 1

    def save_skeleton(self, mesh, colors=None):
        file_name = 'skeleton_%d.obj' % self.skeleton_counter
        if colors:
            write_trimesh_with_color(mesh, os.path.join(self.skeleton_path, file_name), colors)
        else:
            write_trimesh(mesh, os.path.join(self.skeleton_path, file_name))
        self.skeleton_counter += 1

class config:

    config_path = None

    def __init__(self):
        cur_dir = os.path.join(os.path.split(__file__)[0])
        self.config_path = os.path.join(cur_dir, 'config.json')
    
    def load_config(self):
        return json.load(open(self.config_path, 'r'))