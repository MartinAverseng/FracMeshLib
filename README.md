# Generalized meshes


This Matlab toolbox contains implementations for the methods described in the paper <br>
[Averseng, Claeys, Hiptmair: Fractured Meshes. <i>FINEL</i> (2023)](https://www.sciencedirect.com/science/article/pii/S0168874X22001809) <br>
It allows to define and manipulate <i> Generalized meshes</i>, which is a structure modeling meshes of non-manifold geometries. It includes in particular fractured meshes, useful for Finite Element Methods (FEM) in fractured domains, and inflated meshes useful for Boundary Element Methods (BEM) for non-manifold boundaries (see, e.g., [this paper](https://arxiv.org/abs/2310.09204) and [this repository](https://github.com/MartinAverseng/multi-screen-bem3D-ddm/)). It builds on open source toolboxes from [this project](https://github.com/matthieuaussal/gypsilab) for finite element quadratures, assembling and standard mesh management.

| FEM Solution of a 2nd order boundary value problem in a randomly generated "fracture network" with Neumann conditions at the boundary | FEM in space/FD in time solution of the wave equation in the complement a complex obstacle |
:----:|:-----:
|![](doc/fractureNetwork.png)|![](doc/wave.gif)|

Run these examples by opening Matlab in Examples, and run the scripts `fractureNetwork.m` and `wave.m`. Refer to [this code](https://github.com/MartinAverseng/multi-screen-bem3D-ddm/) for BEM computations.
Numerical experiments from the paper about the disk with cracked radius can be reproduced by running the scripts in FEMExperiments folder. 


# Table of contents
1. [What is a Generalized mesh](#what-is-a-generalized-mesh)
2. [Class definition and methods](#class-definition-and-methods)
3. [Finite Element Methods in fractured Meshes](#finite-element-methods-in-fractured-meshes)

## What is a Generalized Mesh

### Definition

A n-dimensional generalized mesh is defined by supplying

- An ordered list <math>L = (S1,S2,...,SN) </math> of n-simplices, allowing repetitions and arbitrary overlaps,
- Some adjacency information between the elements of this list.

The adjacency specifies, for each element i = 1,...,N, and each (n-1)-dimensional face <math>F</math> of the simplex Si, whether i has a neighbor element j <i>through</i> <math>F</math>. It must satisfy the following axioms.

(i) If elements i and j are adjacents through F, then F is a facet of both Si and Sj <br>
(ii) Adjacency is symmetric: element i is adjacent to j through F <i>iff</i> element j is adjacent to i through F<br>
(iii) For each facet <math>F</math> of Si, element i is adjacent to <i>at most</i> one element j through <math>F</math>

### Examples


#### Regular meshes

Every <i>regular</i> mesh (of a manifold geometry) can be regarded as a Generalized mesh, by letting elements i and j be adjacent through F <i>iff</i> S<sub>i</sub> interesected with S<sub>j</sub> is equal to F. Given a regular mesh `m` represented by a Gmsh-like structure (vertices,elements), call the class constructor
```
M = GeneralizedMesh(m)
```
to construct the corresponding generalized mesh.


#### Fractured meshes

Generalized meshes can also represent more complex geometries. <i>Fractured meshes</i>, such as the one represented below are important examples: they correspond to regular n-dimensional meshes in Rn, where some adjacencies have been "dropped" at some (n-1) dimensional interelement interfaces, creating a <i>fracture</i> (red edges below). They can be created by calling
```
M = fracturedMesh(mOmega,mGamma)
```
where `mOmega` is a n-dimensional regular mesh and `mGamma` is a (n-1)-dimensional such whose elements are all faces of some simplex in `mOmega`.


<div>
<img src="doc/FracMeshExample.png" alt="Example of Fractured Mesh" width=500/>
</div>

#### Inflated meshes

Generalized meshes can also represent more exotic geometries, for instance by having several elements on top of each other. An important example is given by so-called <i>inflated meshes</i> which naturally appear as the <i>boundaries</i> of fractured meshes (see more below).


## Class definition and methods

### Properties

The class definition of GeneralizedMesh is

```
classdef GeneralizedMesh

    properties
        vtx; % Nvtx x 3 array of reals (coordinates)
        elt; % Nelt x (n+1) array of numbers in {1,...,Nvtx}
        nei_elt; % Nelt x (n+1) array of numbers in {0,..,Nelt}
        nei_fct; % Nelt x (n+1) array of numbers in {0,...,n+1}
    end

    ...

end
```

The pair (`vtx`,`elt`) works like for a normal Gmsh-type mesh, except now the array elt may have duplicate lines, allowing for distinct elements with the same vertices. The adjacency is encoded via `nei_elt` and `nei_fct`. For a non-zero j,
`nei_elt(i,alpha)=j` and `nei_fct(i,alpha) = beta` means that elements i and j are adjacent through the facet Falpha(i), where Falpha(i) is the (n-1)-subsimplex of Si obtained by removing vertex number alpha. In this case, one has also `nei_elt(j,beta) = i, nei_fct(j,beta) = alpha` and Falpha(i) = Fbeta(j). On the other hand, `nei_elt(i,alpha) = 0` and `nei_fct(i,alpha) = 0` is the convention when element i has no neighbor through Falpha(i).     

### Subsimplices

The d-subsimplices of a generalized mesh `M` are all the d-subsimplices of its elements. Call

```  
[Slist,set2sub,sub2set] = subsimplices(M,d)
```
to get all d-subsimplices of M. The i-th line of Slist, `[Slist(i,1),...,Slist(i,d)]` encodes a unique d-subsimplex S<sub>i</sub> of `M`, referring to vertices by their index in `M.vtx`.

- `set2sub` is an array of size Nelt x Nd, where Nd = (n choose d) is the number of d-subsimplices per element. The k-th line of `set2sub` tells at what positions the subsimplices of element k are in `Slist`.
- `sub2set` is a `Nelt x N` sparse "incidence" matrix filled with 0s and 1s, with coefficient (k,i) equal to 1 when element k contains subsimplex S<sub>i.

### Generalized subfacets

Generalized subfacets are the central concept related to Generalized meshes when it comes to FEM and BEM. For each vertex <math>S</math> in a generalized mesh, one can partition the elements incident to <math>S</math> in connected components, with two elements in the same component if they can be linked by a chain of adjacent elements all incident to <math>S</math>. Each component gamma gives rise to a <i>generalized vertex</i> encoded by a pair <b>s</b> = (S,gamma). The sketch below shows a vertex with 3 associated generalized vertices: <b>s</b><sub>1</sub> = (S,{1}),  <b>s</b><sub>2</sub> = (S,{2,3,4}) and  <b>s</b><sub>3</sub> = (S,{5,6})
<div align="center">
<img src="doc/genVert.jpg" alt="Sketch of generalized vertex" width=600/>
</div>
Generalized d-subfacets are defined analogously. Call

```
[Slist,gamma,I] = M.generalizedSubfacets(d)
```

to compute all generalized d-subfacets of `M`,  {<b>s</b> = (S,gamma)}. Generalized subfacet number i is represented as the pair <b>s</b><sub>i</sub> = (`Slist(i,:)`, `gamma{i}`) where

- `Slist(i,:)`, the i-th line of `Slist`, is the d-subsimplex S
- `gamma{i}` is the list of elements in the component gamma

Moreover, `I(i)` refers to the the d-subsimplex S by its index as returned by `subsimplices(M,d)`.

### Refinement

We provide a method for midpoint refinement of triangular generalized meshes. Refinements for other dimensions have not been implemented yet.
Call

```
M = M.refine(p)
```
to refine `M` p times (thereby reducing the elements diameters by a factor 2<sup>p</sup>)


### Boundary

The boundary of a n-dimensional Generalized Mesh M is a (n-1)-dimensional generalized mesh dM. The definition of the boundary is purely combinatorial, involving again chains of adjacent elements (cf the Fractured Meshes paper for details). If M stands for a normal, regular mesh, then dM stands for the usual boundary of M. The boundary of a generalized mesh is obtained by

```
dM = genBoundary(M)
```
The sketch below illustrates the result when `M` is a fractured mesh as represented on the left. The generalized boundary is represented on the right. Note that in reality, each pair of red edges in the central cross are supposed to be exactly on top of each other (overlayed with the original fracture position), but they have been separated artificially to visualize the adjacency information more clearly. Note that both M and dM have 4 distinct generalized vertices at the cross point.
![](doc/boundary.jpg)

### Intrinsic inflation

When given only the simplices of the boundary of a fractured mesh but no adjacency information (e.g., only the 4 red edges in the left panel of the previous figure), one can in fact recover all adjacency information (i.e., in the previous example, reconstruct the right "inflated" cross) geometrically, via a fast and automatic process. This algorithm is key to perform BEM on non-manifold boundaries dM without need for the mesh M. Given an input Gmsh-like (and possibly non-manifold) mesh `m` representing the geometry of the boundary, either an edge mesh in 2D or a triangular mesh in 3D, call

```
M = intrinsicInflation(m)
```

to use this algorithm. Roughly speaking, `M` represents a two sided version of `m`, with twice as many elements, each corresponding to one orientation of an initial element. The sketch below illustrates the concept.
![](doc/inflation.png)

## Finite Element Methods in fractured Meshes

### Whitney forms on M

Spaces of discrete d-differential forms can be defined on Generalized Meshes, and these can be used to consider conforming Galerkin methods in fractured domains for 2nd order PDEs. The key idea is that the degrees of freedom of the d-dimensional finite element spaces (d= 0: P1 Lagrange, d=1: Nédélec, d=n-1:Raviart-Thomas, d=n: Piecewise constant) are exactly the generalized d-subfacets defined above.

Here we discuss this in the simplest case of d=0. For each generalized vertex <b>s</b> = (S,gamma), the basis function ϕ<sub><b>s</b></sub> is the usual tent function at S <i>multiplied by the indicator function of the union of the elements in gamma</i>. The set Λ<sup>0</sup>(M) of discrete 0-forms, or P1 finite element space, is the vector space spanned by {ϕ<sub><b>s</b></sub>}<sub><b>s</b></sub> for all generalized vertices <b>s</b>.
Given a generalized mesh M, call

```
Lambda0M = GenFem(M,'P1')
```
to construct an object representing Λ<sup>0</sup>(M).


Why this space is important is due to the fact that, when M is a fractured mesh of Omega\\ Gamma, then

Λ<sup>0</sup>(M) = {u in H^1(Omega\\ Gamma) : u is affine on each element of M}.

In other words, Λ<sup>0</sup>(M) is a finite-dimensional subspace of the energy space.



### Assembling of FEM matrices

Consider the PDE

-Δu + cu = f(x) in Ω \\ Γ

with Neumann boundary conditions on ∂Ω and on Γ. The variational formulation of this problem reads

Find u in H<sup>1</sup>( Ω \\ Γ) such that, for all v in H<sup>1</sup>( Ω \\ Γ),<br>
(∇u,∇v)<sub>L<sup>2</sup>(Ω \\ Γ)</sub> + c(u,v)<sub>L<sup>2</sup>(Ω \\ Γ)</sub> = (f,v)<sub>L<sup>2</sup>(Ω \\ Γ)</sub>

Given a fractured mesh M of  Ω \\ Γ, the Galerkin approximation u<sub>h</sub> of the weak solution u is defined by

Find u<sub>h</sub> in Λ<sup>0</sup>(M) such that, for all v<sub>h</sub> in Λ<sup>0</sup>(M),<br>
(∇u<sub>h</sub>,∇v<sub>h</sub>)<sub>L<sup>2</sup>(Ω \\ Γ)</sub> + c(u<sub>h</sub>,v<sub>h</sub>)<sub>L<sup>2</sup>(Ω \\ Γ)</sub> = (f,v<sub>h</sub>)<sub>L<sup>2</sup>(Ω \\ Γ)</sub>

This is equivalent to the linear system

(<b>K</b> + c<b>M</b>)U = L

where the stiffness and mass matrices <b>K</b>, <b>M</b> are given by

<b>K</b><sub>i,j</sub> = (∇ϕ<sub><b>s</b><sub>i</sub></sub>,∇ϕ<sub><b>s</b><sub>j</sub></sub>), &nbsp;&nbsp;&nbsp; <b>M</b><sub>i,j</sub> = (ϕ<sub><b>s</b><sub>i</sub></sub>,ϕ<sub><b>s</b><sub>j</sub></sub>) &nbsp;&nbsp;&nbsp; 1 <= i,j <= N<sub>dof</sub>

where <b>s</b><sub>1</sub>,...,<b>s</b><sub>N<sub>dof</sub></sub> is the set of generalized vertices of M.

To assemble these matrices, first call

```
ngauss = 3;
domOmega =dom(M,ngauss);
```
to construct local Gaussian quadrature rules on each element. Then call
```
Mass = integral(domOmega,Lambda0M,Lambda0M);
K = integral(domOmega,grad(Lambda0M),grad(Lambda0M));
```
where `Lambda0M = GenFem(M,'P1')`. If f(X) = sin(X<sub>1</sub>), then the right hand side vector can be computed via
```
f = @(X)(sin(X(:,1)));
L = integral(domOmega,Lambda0M,f)
```

Solve the linear system using

```
c = 0.01;
U = (K + c*Mass)\L;
```
Plot the solution e.g. using this code
```
[X,T] = Lambda0M.dof;
figure;
patch('Faces',T,'Vertices',X,'FaceVertexCData',U,'FaceColor','interp','EdgeColor','interp');
patch('Faces',mGamma.elt,"Vertices",mGamma.vtx,'LineWidth',3,'EdgeColor','k')
caxis([min(min(U)),max(max(U))]);
colormap(jet);
axis equal
title("Solution u")
hold on

c = colormap;
[c1,c2] = caxis;
nl = size(colormap,1);
l = c1 + (c2 - c1)*(1:nl)/nl;
drawLevelSet(M,U,l,'k');
axis off
title("Plot of solution U")
colorbar;
```
