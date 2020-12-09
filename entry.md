# Libigl-style Tutorial

(Please see [github](https://github.com/estherlin/dual-laplace) for full project and [YouTube]() for video)

While the cotan Laplacian operator has been well studied and implemented for triangle meshes, there is no similar default option for tetrahedral meshes. This is because for triangle meshes, two different approaches from the finite element and finite volume methods both yield the cotan Laplacian. However, for tetrahedral meshes, the finite element method yields the primal Laplacian while the finite volume method yields the dual Laplacian. These two Laplacians are often different and can have different applications. The focus of this project is to construct the dual Laplacian operator for a tetrahedral mesh. 

Starting with a finite volume approach, we consider the action of a Laplacian to be a sparse matrix of weights $w_{ij}$ in:
$$
(\mathbf{L f})_{i}=\sum_{(i, j) \in \mathcal{M}} w_{i j}\left(f_{j}-f_{i}\right) \text{, }\quad w_{i j}=w_{j i}
$$
where $\mathcal{M}$ is our mesh. Applying Stokes theorem to the integrated Laplcian leads to an expression for $w_ij$:
$$
w_{i j}=\frac{\operatorname{Vol}(\star(i, j))}{\operatorname{Vol}(i, j)}
$$
where $\operatorname{Vol}(\cdot)$ is the measure of an element. For example, $\operatorname{Vol}(i,j)$ would be the length of the edge between vertices $i$ and $j$, $\operatorname{Vol}(i,j,k) $ would be area, and $\operatorname{Vol}(i,j,k,l)$ would be volume. Also, $(\star(i,j))$ is the polygon dual to edge $(i, j)$. The vertices of this polygon will be at the circumcentre of the tetrahedron, the two triangles adjacent to edge $(i,j)$, and the midpoint of edge $(i,j)$. Next, we slice the polygon into two triangles. From here, the area of the polygon can be found by determining the signed volume of the tetrahedron formed by each of these triangles with vertex $i$ and dividing it by $\operatorname{Vol}(i,j)$, the height of this tetrahedron. 

We take the approach of summing together the contributions for each triangle, for each tetrahedron in the mesh to construct the dual Laplacian $\mathbf{L}$, similar to how the cotan Laplacian is constructed in libigl. This expression indicates that even for a tetrahedral mesh, the edge weights can be considered on a per triangle basis. The benefit of this approach is that we can construct the diagonal mass matrix simultaneously. The pseudocode used to construct the dual Laplacian $\mathbf{L}$ and the mass matrix $\textbf{M}$ is as follows:

```c++
for (int i: tetrahedrons){
    calculate circumcentre of tet: tet_cc
    for (int j: 4 vertices){
        for (int k: 3 vertices, k != j){
            calculate circumcentre of triangle construct: tri_cc
            calculate midpoint of edge (j,k): edge_cc
            calculate volume of tetrahedron (j, tet_cc, tri_cc, edge_cc)
            calculate weight
              
            L(i,j) += weight;
            L(j,i) += weight;
            L(i,i) -= weight;
            L(j,j) -= weight;
            M(i,i) += volume;
            M(j,j) += volume;
        }
    }
}
```

This code is implemented in ```src/dual_laplacian.h``` and can be called with

```c++
dual_laplacian(TV, TT, L, M);
```

where ```TV``` is our matrix of vertex coordinates and ```TT``` is matrix of tetrahedrons, and L and M are sparse matrix representations of the operators. 

![Solving of a Dirichlet boundary value problem with the dual laplacian](./assets/aggregate_view.png)

​	*Solving of a Dirichlet boundary value problem with the dual laplacian*

![Solving of a Dirichlet boundary value problem with the dual laplacian animation](./assets/dual_cube.gif)

​	*Solving of a Dirichlet boundary value problem with the dual laplacian animation*
