# Simulation Setup
**By Keivan A.Gh.**

## Details
This folder contains all the required files for simulating a case by *rhoReactingCentralFoam* solver. This specific setup is based on the main scenario of this project which is a 70-MPA hydrogen leakage into a semi-confined environment with a downstream obstruction. It is important to note that these are raw files which means no command has been executed in this directory. Therefore, in order to start a simulation you should at least go through these steps:

  **1. `blockMesh`:** Creating the geometry and generating the mesh based on [blockMeshDict](./system/blockMeshDict) file.
  \
  **2. `checkMesh`:** Checking the quality of the generated mesh.
  \
  **3. `setFields`:** Setting the initial conditions of the flow field based on [setFieldsDict](./system/setFieldsDict) file.
  \
  **4. `decomposePar`:** Decomposing the domain into several parts (in this case 16) based on [decomposeParDict](./system/decomposeParDict) file for running a parallel simulation.
  \
  **5. `mpirun -np 16 rhoReactingCentralFoam -parallel`:** Running the simulation by utilizing rhoReactingCentralFoam solver and 16 cores of the CPU.
