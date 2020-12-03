# dual laplace

An implementation of the tetrahedral laplacian operators described in the SGP 2020 paper [Properties of Laplace Operators for Tetrahedral Meshes](https://igl.ethz.ch/projects/LB3D/LB3D.pdf). 

Template copied from [libigl-example-project](https://github.com/libigl/libigl-example-project). 

Author: [esther](http://www.cs.toronto.edu/~lin/)

## Dependencies

The only dependencies are stl, eigen, [libigl](http://libigl.github.io/libigl/) and
the dependencies of the `igl::opengl::glfw::Viewer`.

The cmake build system will attempt to find libigl according to environment variables (e.g., `LIBIGL`) and searching in common desitinations (e.g., `/usr/local/libigl/`). If you haven't installed libigl before, we recommend you to clone a copy of libigl right here:

    cd libigl-example-project/
    git clone https://github.com/libigl/libigl.git

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `example_bin` binary.

## Run

From within the `build` directory just issue:

    ./example

A glfw app should launch displaying a 3D cube.
