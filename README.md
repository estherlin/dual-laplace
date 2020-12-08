# dual laplace

An implementation of the tetrahedral laplacian operators described in the SGP 2020 paper [Properties of Laplace Operators for Tetrahedral Meshes](https://igl.ethz.ch/projects/LB3D/LB3D.pdf). 

In addition to implementing the dyal laplacian, I also provide code written in libigl-style for:
1. ```src/circumcentre3d.h```: calculates the circumcentre of triangles and tetrahedrons in 3D space  
2. ```src/tet_volume.h```: calculates the volume of a tetrahedron

Template from [libigl-example-project](https://github.com/libigl/libigl-example-project). 

Author: [esther]()

## Setup

### Dependencies

-   [libigl](http://libigl.github.io/libigl/)
-   Eigen
-   tetgen

### Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should create a `tetra` binary.

### Run

From within the `build` directory just issue:

    ./tetra
A glfw app should launch displaying a 3D cube.

## Examples

### Laplace equation with Dirchlet boundary conditions

![Cooling one side of a cube](./assets/cube.gif)

In this example, we solve the Laplace equation
$$
\Delta u = 0
$$
where $\Delta$ is the Laplacian operator and $u$ are the values at the vertices of our 3D mesh, subject to Dirichlet boundary conditions
$$
u\vert_{\partial S} = 0

$$

### Smoothing

