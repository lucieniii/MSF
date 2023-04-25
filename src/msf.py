import utils
import trimesh

if __name__ == '__main__':
    logger = utils.logger()
    a = "ss"
    logger.logging(1, a, 'aa', to_console=True)
    mesh = utils.read_mesh_trimesh('models/test/armadillo_21.obj', process=False)
    print(mesh)
    ms = utils.MeshSaver(logger.result_path)
    ms.save_skeleton(mesh)
    ms.save_voronoi(mesh)
    config = utils.config()
    dic = config.load_config()
    print(type(dic), dic)