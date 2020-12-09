# dual laplace

An implementation of the tetrahedral laplacian operators described in the SGP 2020 paper [Properties of Laplace Operators for Tetrahedral Meshes](https://igl.ethz.ch/projects/LB3D/LB3D.pdf). 

In addition to implementing the dual laplacian, I also provide code written in libigl-style for:
1. ```src/circumcentre3d.h```: calculates the circumcentre of triangles and tetrahedrons in 3D space  
2. ```src/tet_volume.h```: calculates the volume of a tetrahedron

The ```.cpp``` files can be found in ```include/```. Documentation on the input/outputs for each function can be found in the header files in ```src/```.

Author: esther, 2020

Project link: [https://github.com/estherlin/dual-laplace](https://github.com/estherlin/dual-laplace)

## Setup

-   Template from [libigl-example-project](https://github.com/libigl/libigl-example-project)
-   Demo adapted from [Example 605 Tetgen](https://github.com/libigl/libigl/blob/master/tutorial/605_Tetgen/main.cpp) and [Example 303 Laplace Equation](https://github.com/libigl/libigl/blob/master/tutorial/303_LaplaceEquation/main.cpp)

### Dependencies

-   [libigl](http://libigl.github.io/libigl/): make sure to have the full install! Have a version of libigl that is built according to the [libigl tutorial](https://libigl.github.io/tutorial/#downloading-libigl)
-   Eigen
-   tetgen: comes with the full install of libigl

### Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should create a `tetra` binary.

### Run

From within the `build` directory, issue:

    ./tetra
A glfw app should launch displaying a 3D cube. 

-   To see cross sections of the cube, enter numbers between 1 and 9
-   To see the entire cube, enter 'B' or 'b'

To see other meshes, issue:

```
./tetra ../data/<file>.off
```

Beware! Other meshes may take a **long** time to load!

## Demo: Laplace equation with Dirchlet boundary conditions

![Cooling one side of a cube](./assets/dual_cube.gif)

In this example, we solve the Laplace equation
$$
\Delta u = 0
$$
where $\Delta$ is the Laplacian operator and $u$ are the values at the vertices of our 3D mesh, subject to Dirichlet boundary conditions
$$
u\vert_{\partial S} = 0
$$



