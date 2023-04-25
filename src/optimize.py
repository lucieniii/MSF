import numpy as np

def get_edges_unique_idx(edges_unique, v1, v2):
    e = np.where((edges_unique[:, 0]==v1)&(edges_unique[:, 1]==v2))
    if e[0].shape[0] == 0:
        e = np.where((edges_unique[:, 0]==v2)&(edges_unique[:, 1]==v1))
    return e[0].item()

def get_adj_face(edges_unique, face_adjacency, v1, v2, face_index):
    e = get_edges_unique_idx(edges_unique, v1, v2)
    faces = face_adjacency[e].flatten()
    return faces[faces!=face_index].item()

def get_adj_face_pair(edges_unique, face_adjacency, v1, v2):
    e = get_edges_unique_idx(edges_unique, v1, v2)
    faces = face_adjacency[e].flatten()
    return faces

def can_be_collapsed(vertex_neighbors, v0, v1):
    n0 = vertex_neighbors[v0]
    n1 = vertex_neighbors[v1]
    return np.intersect1d(n0, n1).shape[0] == 2

def cal_edge_length(vertices, v0, v1):
    return np.linalg.norm(vertices[v0] - vertices[v1])

# v1v0 onto v1v2
def cal_projection(v0, v1, v2):
    v1v2 = v2 - v1
    if np.linalg.norm(v1v2) < 1e-6:
        return v1
    v1v0 = v0 - v1
    projector = v1v2 / np.linalg.norm(v1v2)
    projectee = v1v0
    t = np.dot(projectee, projector)
    return v1 + t * projector

def replace_vertex(face, v_old, v_new):
    new_face = face.copy()
    new_face[face==v_old] = v_new
    return new_face

def cal_edge_length(vertices, v0, v1):
    return np.linalg.norm(vertices[v0] - vertices[v1])

#         v0 ---------- v2
#        /        +     /
#      /      vn       /
#    /     +          /   
#  v1 ---------------v3
def splitter(mesh, D, voronoi_poles, is_fixed, short_edge=0.0):

    v_pairs = [[0,1],[1,2],[0,2]]
    faces = np.asarray(mesh.faces)
    edges_unique = np.asarray(mesh.edges_unique)
    face_adjacency = np.asarray(mesh.face_adjacency)
    face_angles = np.asarray(mesh.face_angles)
    vertices = np.asarray(mesh.vertices)

    f_is_del = np.zeros((mesh.faces.shape[0],), dtype=np.bool_) # f is deleted or not
    f_is_fixed = np.zeros((mesh.faces.shape[0],), dtype=np.bool_)
    f_to_add = np.asarray(mesh.faces) # all v indexs. Note that here it only adds points and dont remove any.
    v_to_add = np.asarray(mesh.vertices) # positions
    vor_to_add = voronoi_poles.copy()
    
    v_count = mesh.vertices.shape[0]
    face_count = faces.shape[0]

    
    for face_index in range(face_count):
        face = faces[face_index]
        if not f_is_fixed[face_index] and max(face_angles[face_index]) > D: # 2.1
            
            # get v0,v1,v2
            v0 = face[np.argmax(face_angles[face_index])] # obtuse_vertex
            v1, v2 = face[face!=v0] # v1,v2 are indexes.  
            
            if is_fixed[v1] and is_fixed[v2]:
                continue
            
            # get indexs of adj faces
            face_adj_index = get_adj_face(edges_unique, face_adjacency, v1, v2, face_index)
            face_adj = faces[face_adj_index]

            if cal_edge_length(vertices, v1, v2) < short_edge:
                continue
            
            if face_angles[face_adj_index][np.where((face_adj!=v1)&(face_adj!=v2))] < D:
                continue              

            # add vertex vn
            v_count += 1
            v_new = cal_projection(*vertices[[v0, v1, v2]])
            v_to_add = np.vstack((v_to_add, v_new))

            vor_new = cal_projection(*vor_to_add[[v0, v1, v2]])
            vor_to_add = np.vstack((vor_to_add, vor_new))

            # delete face
            f_is_del[face_index]= True
            f_is_del[face_adj_index]= True

            # fix adj faces in this iteration
            for i, j in v_pairs:
                f_is_fixed[get_adj_face(edges_unique, face_adjacency, 
                                        face[i], face[j], 
                                        face_index)] = True
                f_is_fixed[get_adj_face(edges_unique, face_adjacency, 
                                        face_adj[i], face_adj[j], 
                                        face_adj_index)] = True

            # # add faces
            vn = v_count - 1
            f_to_add = np.vstack((f_to_add, replace_vertex(face, v1, vn)))
            f_to_add = np.vstack((f_to_add, replace_vertex(face, v2, vn))) # the order matters! if go v0v2vn, not work
            f_to_add = np.vstack((f_to_add, replace_vertex(face_adj, v1, vn))) # why
            f_to_add = np.vstack((f_to_add, replace_vertex(face_adj, v2, vn)))                
            f_is_del = np.hstack((f_is_del, [False,False,False,False]))

    faces_new = f_to_add[~f_is_del]
    v_add = v_count - mesh.vertices.shape[0]

    mesh.vertices = v_to_add
    mesh.faces = faces_new

    assert mesh.is_watertight, 'mesh is not watertight'

    return v_add != 0, v_add, vor_to_add

def collapser(mesh, T, is_fixed):

    faces = np.asarray(mesh.faces.copy())

    edges_length = np.asarray(mesh.edges_unique_length.copy())
    short_idx = edges_length < T
    edges_unique = np.asarray(mesh.edges_unique.copy())[short_idx]
    face_adjacency = np.asarray(mesh.face_adjacency.copy())[short_idx]
    edges_length = edges_length[short_idx]
    vertex_neighbors = mesh.vertex_neighbors.copy()

    ignore = np.zeros((mesh.vertices.shape[0],), dtype=np.bool_)

    f_is_del = np.ones((mesh.faces.shape[0],), dtype=np.bool_)

    c_count = 0

    collapsed = []

    for edge, adj_faces in zip(edges_unique, face_adjacency):
        v0, v1 = edge

        # check can be collapsed
        if ignore[v0] or ignore[v1] or (is_fixed[v0] and is_fixed[v1]) or not can_be_collapsed(vertex_neighbors, v0, v1):
            continue
        if is_fixed[v1]:
            v0, v1 = v1, v0
            
        # reset v0
        mesh.vertices[v0] = 0.5 * (mesh.vertices[v0] + mesh.vertices[v1])

        # delete adj faces
        f0, f1 = adj_faces
        f_is_del[f0] = False
        f_is_del[f1] = False

        # change v1 -> v0
        faces[faces==v1] = v0

        c_count += 1

        collapsed.append(v1)

        ignore[v1] = True
        for nv in vertex_neighbors[v1]:
            ignore[nv] = True

        ignore[v0] = True
        for nv in vertex_neighbors[v0]:
            ignore[nv] = True
        
    mesh.faces = faces[f_is_del]
    mesh.remove_unreferenced_vertices()
    
    return collapsed

def collapse_Edges(mesh, T, voronoi_poles, WL, WH, WM, is_fixed, v_color):
    can_collapse = True
    n_collapsed = 0
    n_it = 0
    while can_collapse:
        collapsed = collapser(mesh, T, is_fixed)
        can_collapse = len(collapsed) != 0
        if can_collapse:
            voronoi_poles = np.delete(voronoi_poles, collapsed, axis=0)
            WL = np.delete(WL, collapsed, axis=0)
            WH = np.delete(WH, collapsed, axis=0)
            WM = np.delete(WM, collapsed, axis=0)
            is_fixed = np.delete(is_fixed, collapsed, axis=0)
            v_color = np.delete(v_color, collapsed, axis=0)
        n_collapsed += len(collapsed)
        n_it += 1
    print(n_collapsed, "vertices collapsed in", n_it, "iterations.")
    assert mesh.is_watertight, 'mesh is not watertight'
    return voronoi_poles, WL, WH, WM, is_fixed, v_color, n_collapsed, n_it

def split_Faces(mesh, D, voronoi_poles, is_fixed, short_edge):
    n_added = 0
    n_it = 0
    can_split = True
    t_is_fixed = is_fixed.copy()
    while can_split:
        can_split, v_count, voronoi_poles = splitter(mesh, D, voronoi_poles, t_is_fixed, short_edge)
        n_added += v_count
        t_is_fixed = np.hstack((is_fixed, np.zeros((n_added,), dtype=bool)))
        n_it += 1
    print(n_added, "vertices added in", n_it, "iterations.")
    assert mesh.is_watertight, 'mesh is not watertight'
    return n_added, voronoi_poles, n_it
