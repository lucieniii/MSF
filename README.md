# UCL-APG-Project
COMP0119: Acquisition and Processing of 3D Geometry Project

A realization of Mean Curvature Skeleton

> **Reference:** Tagliasacchi, A., Alhashim, I., Olson, M. and Zhang, H. (2012), Mean Curvature Skeletons. Computer Graphics Forum, 31: 1735-1744. https://doi.org/10.1111/j.1467-8659.2012.03178.x

## How to use
*If you want, you can adjust the parameters in `src/config.json`.*
```Python
# Assuming you are in the root of project
from src.msf import MSF

# Initial a MSF Instance
msf = MSF()

# Load a mesh
msf.load_mesh("../models/armadillo.obj")

# Let's iterate! Of course you should repeate it.
msf.iterate()

# Then see the skeletons in the latest folder of the `result` folder.
```

## Parameter Explaination
*Set the parameters in `src/config.json`*

+ `wL`, wL in the paper;
+ `wH`, wL in the paper;
+ `wM`, wL in the paper;
+ `use_dynamic_scale`, If `true`, then allow adjusting `scale` between iterations. The `scale` will change between `scale_min` and `scale_max` repeatedly at intervals of `scale_delta`, which can speed up convergence;
+ `scale`, The $\varepsilon$ in the paper will be scale multiply the diagonal length of the mesh's bounding box;
+ `scale_min`, Only be used when `use_dynamic_scale=true`;
+ `scale_delta`, Only be used when `use_dynamic_scale=true`;
+ `scale_max`, Only be used when `use_dynamic_scale=true`;
+ `scale_fix`, Be used to fix vertices;
+ `alpha`, The $\theta$ in the paper;
+ `smooth_lam`, $\lambda$ in the Implict Laplacian Smoothing (step size);
+ `smooth_it`, The iteration number of the Implict Laplacian Smoothing;
+ `use_reconstruction_in_voronoi`, If `true`, use mesh reconstruction before calculating Voronoi poles, which can produce a better approximation of medial axis;
+ `k`, Parameter in mesh reconstruction, only be used when `use_reconstruction_in_voronoi=true`, which means using the eigen vectors corresponding to the `k` smallest eigen values to perform reconstruction;
+ `Laplacian_type`, Determine which Laplace–Beltrami operator to use. Choose one of: `uniform`, `cotangent` and `tangent`.
+ `lsqr_args`, Parameters for `scipy.sparse.linalg.lsqr`. Generally not modified.

## Algorithm Procedure Code
![procedure_code](docs/procedure_code.png)