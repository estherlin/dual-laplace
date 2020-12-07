# dual laplace

An implementation of the tetrahedral laplacian operators described in the SGP 2020 paper [Properties of Laplace Operators for Tetrahedral Meshes](https://igl.ethz.ch/projects/LB3D/LB3D.pdf). 

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

