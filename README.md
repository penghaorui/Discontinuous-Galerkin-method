# Discontinuous-Galerkin-method
Implementation of the Discontinuous Galerkin method for elastic modeling wave propagation in complex media with triangluar meshes using MPI and OpenMP.  

Author information: **Haorui Peng**  
This code was created during my master project in China University of Geosciences(Wuhan), from 2014 to 2017.  

The implementation is mainly based on the paper below, but modified with Central Flux.   

XUE Zhao, DONG Liang-Guo, LI Xiao-Bo, LIU Yu-Zhu. Discontinuous Galerkin finite-element method for elastic wave modeling including surface topography, Chinese Journal of Geophysics (in Chinese), 2014, 57(4): 1209-1223, doi: 10.6038/cjg20140418.  

The application example of this code is as below:   
Haorui Peng et al., Preliminary application of discontinuous Galerkin method in seismic response of fractured-vuggy 
carbonatite reservoirs, SPG/SEG International Geophysical Conference, 2016. 

To generate meshes for a given model, please check  
https://www.cs.cmu.edu/~quake/triangle.html 

Here is an example of the triangular meshing for a slope model with coupled acoustic-elastic interface:\
<img src="https://github.com/penghaorui/Discontinuous-Galerkin-method/blob/main/input_files/mesh/slope_model_mesh.png" width="600" height="600" />

Two versions are provided here:
 > - a) using MPI and OpenMP
 > - b) using only OpenMP

The input_files directory provides the files for the input. Two subfolders are included in the input_file directory. The input_files/mesh contains the mesh files and the matlab codes to generate them.

input_files/precomputed_matrix/ contains the matlab codes to generate precomputed matrices that will be used as input to the DG code later. The precomputed matrices are necessary according the DG theory. Part of these matlab codes are given by the book: Nodal Discontinuous Galerkin Methods: Algorithms, Analysis, and Applications.

The output_files contain the files for the output. \
Here is an example of the output wavefield snapshot for the slope model
<img src="https://github.com/penghaorui/Discontinuous-Galerkin-method/blob/main/output_files/slope_Vx_snapshot.jpeg" width="600" height="600" />

To use the code follow the steps below:  
step 1: create mesh files by running input_files/mesh/**main.m**  
step 2:  
>if choose to use MPI with 2 cores, run the command below:
> - mpicc -o dg DG_leap_Frog_MPI.c -lm -fopenmp
> - mpirun -np 2 dg

>if choose only to use OpenMP, run the command below:
> - gcc -o dg DG_leap_Frog_OpenMP.c -lm -fopenmp
> - ./dg
