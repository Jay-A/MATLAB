"""""""""""""""""
LocalTFEM (incomplete README...)
"""""""""""""""""
...................................................................................................
A continuous 2D finite element for unstructured, simplicial meshes requiring only local operations
...................................................................................................

.. contents:: Overview
   :depth: 3

===================
About
===================

This is a *(p-1)*-exact, continuous, triangular finite element method for explicit time advancement
This version is meant to project 2nd order polynomials exactly and achieve the expected 3rd order convergence in space.
The project is based on article currently availabe as preprint at https://arxiv.org/abs/1906.10774

===================
Dependencies
===================

There aren't really any dependencies. The project is written using MATLAB 2020a and makes use of 
symbolic variables to define the modified Dubiner mass matrix.

===================
Usage
===================

The project is structured as follows

::

    MATLAB-LocalTFEM/
    │
    ├── README.rst
    │ 
    ├── test_p345Conv.m
    ├── test_DubinerSquareHeatedPlate.m
    │
    ├── private/          
    │   ├── assembly/
    │   │   ├── Create_Global_Addresses.m
    │   │   ├── Create_Global_Mass.m
    │   │   ├── Create_Global_RHS.m
    │   │   ├── Create_Global_RHSLT.m
    │   │   └── Create_Map.m
    │   ├── dubiner/
    │   │   ├── Create_Poly.m
    │   │   ├── Dubiner_Mass.m
    │   │   ├── Dubiner_Modes.m
    │   │   ├── GradPhi.m
    │   │   └── Jacobi_Coeffs.m
    │   ├── libraries/
    │   │   └── ...
    │   ├── local/
    │   │   ├── Eval_RHS.m
    │   │   └── Eval_RHSLT.m
    │   ├── mesh/
    │   │   ├── mshLoader.m
    │   │   └── structured_mesh.m
    │   ├── quadrature/
    │   │   ├── lglnodes.m
    │   │   ├── lgrnodes.m
    │   │   ├── lgwt.m
    │   │   └── RefTri_Quad.m
    │   └── visualization/
    │       ├── Make_Video.m
    │       └── View_Soln.m
    └──


Currently test_p345Conv.m evaluates the L2 errors of the projection of F0 (line 7) onto the 
locally defined basis for p=3, 4, or 5 (line 10) over a sequence of structured meshes ranging 
in size from m_first to m_final (lines 25 and 26). Running test_p345Conv.m returns 
Mesh information ( Areas, Boundaries, Edges, Elements, Nodes, Verts, map ), Basic Information (p, F0),
approximate solution (ULT), and the errors in a matrix:
 
L2_Errors.mat a table of L2 errors for given mesh sizes

======  ======    ======
 999    L2 Err     999
------  ------    ------
(step)             h 
======  ======    ======
 0       err0      h0
 1       err1      h1
 2       err2      h2
 ...     err3      ...
======  ======    ======

----------------------
Tests
----------------------

Current tests available:

Convergence of LT for p=3,4,5
------------------------------
The MATLAB script test_p345Conv.m creates convergence plots to demonstrate the desired convergence rates of the novel lower-triangular pseudo-mass matrix method.
This method is shown to be *(p-1)*-exact and achieves the *p-th* order spatial convergence rates.

Heat Equation
----------------------
The MATLAB script test_DubinerSquareHeatedPlate.m creates an .avi video demonstrating a heated plate with one homogeneous Dirichlet boundary condition on the left.
The script allows the user to change initial heat distribution, mesh size, the polynomial degree of the modified Dubiner basis used, final time, time step size, and frame quality for the video frame rendering.

--------------------------------
Local Integration for LT scheme
--------------------------------
Not yet.


===================
Acknowledgements
===================

Probably Helenbrook
