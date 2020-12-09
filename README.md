# dual laplace

An implementation of the tetrahedral laplacian operators described in the SGP 2020 paper [Properties of Laplace Operators for Tetrahedral Meshes](https://igl.ethz.ch/projects/LB3D/LB3D.pdf). 

In addition to implementing the dual laplacian in ```src/dual_laplacian.h```, I also provide code written in libigl-style for:
1. ```src/circumcentre3d.h```: calculates the circumcentre of triangles and tetrahedrons in 3D space  
2. ```src/tet_volume.h```: calculates the volume of a tetrahedron

Maybe these can also be useful additions to libigl!

Please refer to:

-   ```entry.md``` for the libigl tutorial-style explanation of this project (also submitted on Markus)
-   [YouTube](https://youtu.be/0IMYSC3cPkQ) for a 1 min demo of this project
-   [GitHub (this project)](https://github.com/estherlin/dual-laplace) for everything else

## Structure of Repo

-   ```include/```: contains header files with documentation
-   ```src/```: contains ```.cpp``` files
-   ```data/```: contains ```.off``` files used in demo and development
-   ```assets/```: contains screenshots and ```.gif``` files demo-ing the program

Links to mathematical equations and derivations are cited in the header and/or ```.cpp``` files. 

## Setup

**Credits to**:

-   Template taken from [libigl-example-project](https://github.com/libigl/libigl-example-project): this was a great place to start!
-   Demo adapted from [Example 605 Tetgen](https://github.com/libigl/libigl/blob/master/tutorial/605_Tetgen/main.cpp) and [Example 303 Laplace Equation](https://github.com/libigl/libigl/blob/master/tutorial/303_LaplaceEquation/main.cpp)
-   The header file for the dual laplacian was taken from [assignment Github](https://github.com/alecjacobson/geometry-processing-libigl-implementations/blob/main/Laplace%20Operators%20for%20Tetrahedral%20Meshes/dual_laplacian.h)

### Dependencies

-   [libigl](http://libigl.github.io/libigl/): make sure to have the full install! Have a version of libigl that is built according to the [libigl tutorial](https://libigl.github.io/tutorial/#downloading-libigl)
-   [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
-   [Tetgen](http://wias-berlin.de/software/tetgen/): comes with the full install of libigl. This can be a ***massive*** pain to install on its own, it's much easier to just go through the entire libigl installation.

### Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should create a `tetra` binary.

### Run

From within the `build` directory, run with:

    ./tetra
A glfw app should launch displaying a 3D cube. 

-   To see cross sections of the cube, enter numbers between 1 and 9
-   To see the entire cube, enter 'B' or 'b'

To see other meshes, issue:

```
./tetra ../data/<file>.off
```

Beware! Other meshes may take a ***long*** time to load!

## Background

While the cotan Laplacian operator has been well studied and implemented for triangle meshes, there is no similar default option for tetrahedral meshes. This is because for triangle meshes, two different approaches from the finite element and finite volume methods both yield the cotan Laplacian. However, for tetrahedral meshes, the finite element method yields the primal Laplacian while the finite volume method yields the dual Laplacian. These two Laplacians are often different and can have different applications. The focus of this project is to construct the dual Laplacian operator for a tetrahedral mesh. 

Starting with a finite volume approach, we consider the action of a Laplacian to be a sparse matrix of weights $w_{ij}$ in:
$$
(\mathbf{L f})_{i}=\sum_{(i, j) \in \mathcal{M}} w_{i j}\left(f_{j}-f_{i}\right) \text{, }\quad w_{i j}=w_{j i}
$$
where $\mathcal{M}$ is our mesh. Applying Stokes theorem to the integrated Laplacian leads to an expression for $w_ij$:
$$
w_{i j}=\frac{\operatorname{Vol}(\star(i, j))}{\operatorname{Vol}(i, j)}
$$
where $\operatorname{Vol}(\cdot)$ is the measure of an element. For example, $\operatorname{Vol}(i,j)$ would be the length of the edge between vertices $i$ and $j$, $\operatorname{Vol}(i,j,k) $ would be area of a triangle, and $\operatorname{Vol}(i,j,k,l)$ would be volume of a tetrahedron. Also, $(\star(i,j))$ is the polygon dual to edge $(i, j)$. The vertices of this polygon will be at (1) the circumcentre of the tetrahedron: $\textbf{c}_{ijkl}$, (2-3) the circumcentre of the two triangles adjacent to edge $ (i,j)$: $\textbf{c}_{ijk}, \textbf{c}_{ijl}$, and (4) the midpoint of edge $(i,j)$: $\textbf{c}_{ij}$. Next, we slice the polygon into two triangles. From here, the area of the polygon can be found by determining the signed volume of the tetrahedron formed by each of these triangles with vertex $i$ and dividing it by $\frac{1}{2}\operatorname{Vol}(i,j)$, the height of this tetrahedron. Finally, we get an expression for $w_{ij}$ in terms of vertices:
$$
w_{i j}=\sum_{(i, j, k, l)} w_{i j k l}+w_{i j l k}
$$
where
$$
w_{i j k l}=6 \frac{\operatorname{Vol}\left(j,\mathbf{c}_{i j}, \mathbf{c}_{i j k}, \mathbf{c}_{i j k l}\right)}{\operatorname{Vol}(i, j)^{2}}\text{ , }\quad w_{i j lk}=6 \frac{\operatorname{Vol}\left(j,\mathbf{c}_{i j}, \mathbf{c}_{i j kl}, \mathbf{c}_{i j  l}\right)}{\operatorname{Vol}(i, j)^{2}}
$$
This expression for the weights indicates that even for a tetrahedral mesh, the edge weights can be considered on a per triangle basis. We take the approach of summing together the contributions for each triangle, for each tetrahedron in the mesh to construct the dual Laplacian $\mathbf{L}$, similar to how the cotan Laplacian is constructed in my smoothing assignment. The benefit of this approach is that we can construct the diagonal mass matrix simultaneously, where the entries are the volume contributions $\operatorname{Vol}\left(j, \mathbf{c}_{i j}, \mathbf{c}_{i j k}, \mathbf{c}_{i j k l}\right)$.

### Implementation

The pseudocode used to construct the dual Laplacian $\mathbf{L}$ and the mass matrix $\textbf{M}$ is as follows:

```c++
for (int i: tetrahedrons){
    calculate circumcentre of tet: tet_cc
    for (int j: 4 vertices){
        for (int k: 3 vertices, k != j){
            calculate circumcentre of triangle construct: tri_cc
            calculate midpoint of edge (j,k): edge_cc
            calculate volume of tetrahedron (j, tet_cc, tri_cc, edge_cc): vol
            calculate weight: w
              
            L(j,k) += w;
            L(k,j) += w;
            L(j,j) -= w;
            L(k,k) -= w;
            M(j,j) += vol;
            M(k,k) += vol;
        }
    }
}
```

This code is implemented in ```src/dual_laplacian.h``` and can be called with

```c++
dual_laplacian(TV, TT, L, M);
```

where ```TV``` is our matrix of vertex coordinates and ```TT``` is matrix of tetrahedrons, and L and M are sparse matrix representations of the operators. 

To help with calculations, I implemented two addition functions in libigl-style:

1.  ```src/circumcentre3d.h```: calculates the circumcentre of triangles and tetrahedrons in 3D space. In the process of implementing this, I found an excellent discussion on all things circumcentre from [Jonathan Shewchuk](https://people.eecs.berkeley.edu/~jrs/) [here](https://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html).  Through his work (and my debugging :D), I discovered that while solutions with cross products are neat, cross products norms in the denominator can be unstable. This manifested in overflow problems for me. 
2.  ```src/tet_volume.h```: calculates the volume of a tetrahedron

## Demo

### Overview

In this example, we solve the Laplace equation on a unit cube:
$$
\Delta u = 0
$$
where $\Delta$ is the Laplacian operator and $u$ are the values at the vertices of our 3D mesh, subject to Dirichlet boundary conditions
$$
u\vert_{\partial S} = 0
$$

where $\partial S$ is the purple face of the cube. The example is adapted from libigl [Example 303 Laplace Equation](https://github.com/libigl/libigl/blob/master/tutorial/303_LaplaceEquation/main.cpp) and the cross sectional views are adapted from libigl [Example 605 Tetgen](https://github.com/libigl/libigl/blob/master/tutorial/605_Tetgen/main.cpp). 

![Solving of a Dirichlet boundary value problem with the dual laplacian](assets/aggregate_view.png)

​	*Solving of a Dirichlet boundary value problem with the dual laplacian on the unit cube. A face is held at zero (purple).*

### Results

In the following GIF, I show the results of holding a face of a unit cube at Dirichlet boundary conditions (set to 0). The cross sections show the  gradient of values changing over the vertices. 

![Cooling one side of a cube](./assets/dual_cube.gif)

​	*Solving of a Dirichlet boundary value problem with the dual laplacian animation*

I'm quite suspicious of these results. While I have confidence in my implementation of the dual Laplacian (see testing section below), the results do not make physical sense for the following reasons:

-   There seems to be a discontinuity/rapid change close to the face held at zero and the face opposite. Near the face held at zero, I expect the interior values to also be close to zero (purplish), but they are not. Similarly, near the face on the opposite end, I expect the interior values to be close to yellow, but they are not. 
-   Very few boundary effects. This is especially apparent in the GIF below. For such a coarsely spaced mesh, I was expecting spikes and aberrations at the corners. While the corners and edges are different from the interior, the difference seems too small. 

![Solving of a Dirichlet boundary value problem with the dual laplacian overhead animation](./assets/dual_overhead_view.gif)

​	*Solving of a Dirichlet boundary value problem with the dual laplacian overhead animation*

-   If we solve the same problem, but using the primal Laplacian, we see similar results with the last two points apparent as well.

    ![Solving of a Dirichlet boundary value problem with the primal laplacian](./assets/primal_cube.gif)

    *Solving of a Dirichlet boundary value problem with the primal Laplacian.* 

These physical inconsistencies may be attributed to a coarse meshing, but I think there is some underlying mistake in the linear systems solving that I neglected. A finer mesh can be used to investigate this, but unfortunately I ran out of time (both for computations and project management :( ) The coarse mesh is also responsible for that one missing tetrahedron in the last face. 

### Testing

The authors of the paper already have [code](https://igl.ethz.ch/projects/LB3D/dualLaplace.cpp) on their [project website](https://igl.ethz.ch/projects/LB3D/). I viewed their code for validation purposes. I ran their code on my unit cube and compared my $\mathbf{L}$, $\mathbf{M}$ operators with theirs, as well as the calculations in ```src/tet_volume.cpp``` and ```src/circumcentre3d.cpp``` with theirs. While I think their code is more elegant, mine does what it needs to. :) Also, the helper functions work too!

## Final Comments

This was a very fun (and contained) paper to implement! A couple of things I would have changed:

-   Implemented a demo that utilizes the mass matrix that was also built (in an application like smoothing)
-   Experimented with more shapes. I found this particularly difficult because meshing time can be dependent on the geometry and volume of the shape. 
-   Learned more about creating meshes that are good for tetrahedralizing. Most of the time spent on this project was on struggles with getting Tetgen set up, Tetgen runtimes, converting between file types, and wrangling non-Delaunay meshes. 

Thank you [CSC2530 Team](https://github.com/alecjacobson/geometry-processing-csc2520) for delivering such a fun, applicable, and hands-on course!









