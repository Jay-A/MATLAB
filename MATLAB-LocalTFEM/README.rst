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

This is a $(p-1)$-exact, continuous, triangular finite element method for explicit time advancement
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
    │   │   └── lglnodes.m
    │   └── visualization/
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

All tests are currently in root directory

----------------------
Heat Equation
----------------------

Not yet

Visualizations
--------------------------

Not yet

===================
Acknowledgements
===================

Probably Helenbrook
