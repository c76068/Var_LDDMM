# Var_LDDMM

This is a MATLAB implementation of the papers [Diffeomorphic Registration of Discrete Geometric Distributions](https://www.worldscientific.com/doi/abs/10.1142/9789811200137_0003) and [Metrics, quantization and registration in varifold spaces](https://arxiv.org/abs/1903.11196) by Hsi-Wei Hsieh and [Nicolas Charon](http://www.cis.jhu.edu/~charon/). The package provides tools to permorm diffeomorphic registration, interpolation and compression of geometric shapes (such as point clouds, discrete cruves or triangulated surfaces) that are represented as ***discrete varifolds***.

## References
If you use this code for your research, please cite our papers:

```
@article{hsieh2019diffeomorphic,
  title={Diffeomorphic Registration of Discrete Geometric Distributions},
  author={Hsieh, Hsi-Wei and Charon, Nicolas},
  journal={Mathematics Of Shapes And Applications},
  volume={37},
  pages={45},
  year={2019},
  publisher={World Scientific}
}

@article{hsieh2019metrics,
  title={Metrics, quantization and registration in varifold spaces},
  author={Hsieh, Hsi-Wei and Charon, Nicolas},
  journal={arXiv},
  pages={arXiv--1903},
  year={2019}
}
```

## Set up
To run functions without GPU acceleration, you will only need MATLAB and the [HANSO](https://cs.nyu.edu/overton/software/hanso/) library for L_BFGS optimization. We have included the HANSO routines in the folder `/src/optimization/` for convenience.

Basic dependencies:
* MATLAB
* [HANSO version 2.2](https://cs.nyu.edu/overton/software/hanso/): in `/src/optimization/`

Var_LDDMM also includes GPU acceleration implementations using the [KeOps](https://www.kernel-operations.io/keops/matlab/index.html) library (libkeops-master in folder '/src/'). It must be used on a machine equipped with an NVIDIA graphics card with recent CUDA drivers installed.

Optional dependencies (GPU acceleration):
* [KeOps](https://www.kernel-operations.io/keops/matlab/index.html): libkeops-master in folder '/src/' 

## Usage
### Registration:
The approach relies on the Large Deformation Diffeomorphic Metric Mapping (LDDMM) model with a *geodesic shooting* scheme. The cost function, which is the sum of a deformation penalty and a kernel fidelity metric between the transformed source and the target varifolds, is minimized with respect to the ***momenta*** variables of the deformation.  
Use `registration` function for LDDMM varifolds matching algorithm:
```Matlab
[P_op,summary]= registration(Source,Target,defo,objfun,options)
```
- Input:
  - `Source`/`Target`: source/target shape represented as a discrete varifold (structure):
    - `Source.center`: a `N-by-m` array contains positions of `N` diracs in dimension `m` 
    - `Source.vector`: a cell array with length `d` carries a frame for the grassmanian, each entry `Source.vector{i}` is a `N-by-m` array
  - `defo`: options for deformation and numerical ODE:
    - `defo.kernel_size_mom`: deformation kernel size
    - `defo.nb_euler_steps`: number of time steps in the forward and backward integration
    - `defo.odemethod`: numerical scheme for ODE integration, possible values: `'rk4'` or `'middle_point'`
    - `defo.method`: compute the deformation kernel operations using only matlab or GPU with keops, possible values: `'keops'` or `'matlab'`

  - `objfun`: options for varifold data attachment term:
    - `objfun.kernel_geom`: spatial kernel type, possible values:`'gaussian'` or `'cauchy'`
    - `objfun.kernel_size_geom`: spatial kernel sizes in an array `[a_1,...,a_K]`, successively runs K times optimizations with kernel sizes `a_1,...,a_K`, with each optimization using the estimated momemta from the previous run as initialization. The length of the array should be consistent with `kernel_size_grass`
    - `objfun.kernel_grass`: Grassmanian/orientation kernel type, possible values: `'linear'`, `'gaussian_oriented'`, `'gaussian_unoriented'`, `'binet'`
    - `objfun.kernel_size_grass`: Grassmanian/orientation kernel sizes in an array `[b_1,...,b_K]`
    - `objfun.method`='keops': compute kernel operation in data attachment term using only matlab or GPU with keops, possible values: `'keops'` or `'matlab'`
    - `objfun.lambda`: weight parameter in front of the data attachment term

  - `options`: options for L-BFGS:
    - `options.record_history`: `true` or `false`
    - `options.maxit`: maximum number of iterations(default `1000`) (applies to each starting vector)
    - `options.nvec`: `0` for full BFGS matrix update, otherwise specifies number of vectors to save and use in the limited memory updates (default: `0`)
    - `options.prtlevel`: one of `0` (no printing), `1` (minimal), `2` (verbose) (default: `1`)

- Output:
  - `P_op`: the optimized momenta stored in a structure
  - `summary`: summary of the optimization procedure

### Compression/Quantization:
Use the function `proj2_M_dirac` to compress a discrete varifold to a sparser discrete varifold with exactly `M` diracs:  
```Matlab
[Y,summary] = proj2_M_dirac(X,X_ini,objfun,options)
```
- Input:
  - `X`: full discrete varifold to be compressed
  - `X_ini`: initial discrete varifold with `M` diracs
  - `objfun` and `options`: same as in `registration`, the field `objfun.lambda` is not needed in compression.
- Output:
  - `Y`: compressed varifold
  - `summary`: summary of the optimization procedure
  
### Examples:
See the two script files in the `Demo scripts` folder for some examples of basic use of the code. The first script `script_Bone_Bottle_quantization_registration.m` computes the compressions and registrations of two curves for M= 25, 40, and 150.
<table align='center'>
<tr align='center'>
<td> Relative quantization errors </td>
<td> M=25</td>
<td> M=40</td>
<td> M=150</td>
</tr>
<tr>
<td><img src = 'images/Bone_Bottle_quantization_energy.png'>
<td><img src = 'images/Bottle_quan25.png'>
<td><img src = 'images/Bottle_quan40.png'>
<td><img src = 'images/Bottle_quan150.png'>
</tr>
<tr>
<td> Difference to thegroundtruth optimal energy </td>
<td> M=25</td>
<td> M=40</td>
<td> M=150</td>
</tr>
<tr>
<td><img src = 'images/Bone_Bottle_quantization_energy.png'>
<td><img src = 'images/Bottle_quan_matching25.png'>
<td><img src = 'images/Bottle_quan_matching40.png'>
<td><img src = 'images/Bottle_quan_matching150.png'>
</tr>
</table>

The second script `script_amygdala.m` registers two amygdala surfaces using GPU acceleration:
<table align='center'>
<tr align='center'>
<td> Amygdala surfaces</td>
<td> Varifold representation</td>
</tr>
<tr>
<td><img src = 'images/amyg_surface.gif'>
<td><img src = 'images/amyg_varifold.gif'>
</tr>
</table>

## Licence

This program is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.

This program is distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with 
this program. If not, see http://www.gnu.org/licenses/.


## Contacts

* Nicolas Charon (charon@cis.jhu.edu)
* Hsi-Wei Hsieh (hhsieh9@jhu.edu)
