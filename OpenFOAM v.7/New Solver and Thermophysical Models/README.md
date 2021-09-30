# Solver and Thermophysical Models
**By Keivan A.Gh.**

This folder contains all the required OpenFOAM utilities that have been used in [this master's thesis](https://github.com/Keivan-A-Gh/Master_Thesis/blob/main/Documents/Master's%20Thesis%20-%20Final%20Edition.pdf).


## Main Files
### 1. rhoReactingCentralFoam Solver
The [utilized solver](https://github.com/Keivan-A-Gh/Master_Thesis/tree/main/OpenFOAM%20v.7/New%20Solver%20and%20Thermophysical%20Models/applications/solvers/compressible/rhoReactingCentralFoam) in this thesis is generally based on [rhoReactingCentralFoam](https://github.com/duncanam/thesis) solver written by Duncan A. McGough at the University of Colorado Boulder, although there are some minor changes and additions in this version to make the original solver compatible with our needs and also the version of the utilized OpenFOAM.

### 2. New Thermodynamic and Transport Models
In this project, for calculating thermodynamic properties (C<sub>p</sub>,  H, S, and etc.) a new thermodynamic model called [nasa9poly](https://github.com/Keivan-A-Gh/Master_Thesis/tree/main/OpenFOAM%20v.7/New%20Solver%20and%20Thermophysical%20Models/src/thermophysicalModels/specie/thermo/nasa9Poly) has been implemented which is based on [NASA-9](https://github.com/Keivan-A-Gh/Master_Thesis/blob/main/OpenFOAM%20v.7/New%20Thermophysical%20Coefficients/2002_(NASA%209)_McBride%2C%20Gordon_NASA%20Glenn%20Coefficients%20for%20Calculating%20Thermodynamic%20Properties%20of%20Individual%20Species.pdf) polynomials for individual species. Furthermore, for calculating transport properties (κ, μ, and etc.) a new transport model called [nasaPoly](https://github.com/Keivan-A-Gh/Master_Thesis/tree/main/OpenFOAM%20v.7/New%20Solver%20and%20Thermophysical%20Models/src/thermophysicalModels/specie/transport/nasaPoly) has been implemented which is based on [NASA-7](https://github.com/Keivan-A-Gh/Master_Thesis/blob/main/OpenFOAM%20v.7/New%20Thermophysical%20Coefficients/1993_(NASA%207)_McBride%2C%20Gordon_Coefficients%20for%20Calculating%20Thermodynamic%20and%20Transport%20Properties%20of%20Individual%20Species.pdf) transport polynomial for individual species.

### 3. New Thermophysical Models
As discussed in the [report](https://github.com/Keivan-A-Gh/Master_Thesis/blob/main/Documents/Master's%20Thesis%20-%20Final%20Edition.pdf), two new thermophysical models have been implemented for this project to enable OpenFOAM to use a non-ideal-gas equation of state along with the new solver as well as the new thermodynamic and transport models to simulate the main scenario. All the required entries that have been added to the specific paths in *basic*, *reactionThermo*, and *specie* libraries are explained in *Chapter 4.3* of the [report](https://github.com/Keivan-A-Gh/Master_Thesis/blob/main/Documents/Master's%20Thesis%20-%20Final%20Edition.pdf). Generally, these two thermophysical models are as follows:

Entries |Initial Model|Final Model
| :--- | :---: | :---:
**`type`**            | `hePsiThermo`            | `hePsiThermo` 
**`mixture`**         | `reactingMixture`        | `reactingMixture` 
**`transport`**       | `sutherland`             | `nasaPoly` 
**`thermo`**          | `janaf`                  | `nasa9Poly` 
**`equationOfState`** | `PengRobinsonGas`        | `PengRobinsonGas` 
**`specie`**          | `specie`                 | `specie` 
**`energy`**          | `sensibleInternalEnergy` | `sensibleInternalEnergy`


## Manual
In order to use the solver, you need to download the [applications](https://github.com/Keivan-A-Gh/Master_Thesis/tree/main/OpenFOAM%20v.7/New%20Solver%20and%20Thermophysical%20Models/applications) folder and merge it to the *applications* folder in your OpenFOAM user working directory (your working directory is located under this path: ``` $WM_PROJECT_USER_DIR ```). After that, you can simply go to the [rhoReactingCentralFoam](https://github.com/Keivan-A-Gh/Master_Thesis/tree/main/OpenFOAM%20v.7/New%20Solver%20and%20Thermophysical%20Models/applications/solvers/compressible/rhoReactingCentralFoam) folder and compile the solver by using the `wmake` command.

Also, for using the new thermodynamic and transport models as well as the new thermophysical models, you just need to download the [src](https://github.com/Keivan-A-Gh/Master_Thesis/tree/main/OpenFOAM%20v.7/New%20Solver%20and%20Thermophysical%20Models/src) folder and as before move it to your OpenFOAM user working directory. This time, first you need to compile the [basic](https://github.com/Keivan-A-Gh/Master_Thesis/tree/main/OpenFOAM%20v.7/New%20Solver%20and%20Thermophysical%20Models/src/thermophysicalModels/basic) library, then the [reactionThermo](https://github.com/Keivan-A-Gh/Master_Thesis/tree/main/OpenFOAM%20v.7/New%20Solver%20and%20Thermophysical%20Models/src/thermophysicalModels/reactionThermo) and [specie](https://github.com/Keivan-A-Gh/Master_Thesis/tree/main/OpenFOAM%20v.7/New%20Solver%20and%20Thermophysical%20Models/src/thermophysicalModels/specie) libraries; failure to follow this sequence will lead to an error in the compilation.


## Citation
For citation, please refer to the [README](https://github.com/Keivan-A-Gh/Master_Thesis/blob/main/README.md) file in the main folder.
