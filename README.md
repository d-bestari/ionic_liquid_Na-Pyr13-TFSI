# MOLECULAR DYNAMICS SIMULATION OF THE EFFECT OF NaTFSI CONCENTRATION AND TEMPERATURE ON IONIC LIQUID Pyr13TFSI 
This repository contains molecular dynamics (MD) simulation data and analysis 
for ionic liquids Na-TFSI and Pyr13-TFSI with varying salt concentrations and temperatures.
The simulations were performed using **GROMACS** with CL&Pol force field parameters.

⚠️ This project is under development. Features may change at any time.


## Computational Resources  
Simulations were performed on the [High Performance Computing (HPC) at BRIN](https://hpc.brin.go.id/) (Badan Riset dan Inovasi Nasional), Indonesia.

## Requirements

### Materials / Force Field
- [CL&P Force Field](https://github.com/paduagroup/clandp) – Force field for ionic liquids from Padua Group

### Simulation Software
- [GROMACS](https://www.gromacs.org) – Molecular dynamics simulation software
- [FFTool](https://github.com/paduagroup/fftool) – Tool for generating initial configurations and force field parameters
- [Packmol](http://www.ime.unicamp.br/~martinez/packmol/) – Tool for building initial molecular boxes

### Scripting / Analysis
- [Python](https://www.python.org) – Programming language
  - NumPy – Numerical computations


## Simulation Overview
- **Ionic liquids**: Na-TFSI and Pyr13-TFSI  
- **Salt concentration variations**: 0.0, 0.1, 0.2, 0.3, 0.4 mole fraction  
- **Temperatures**: 300 K, 350 K, 400 K

### Simulation Workflow
1. **Box Volume Determination**
   - Energy Minimization → Annealing → NVT Equilibration → NPT Equilibration → NPT Equilibration cont → NPT Production

3. **Production NVT**
   - Uses volume obtained from NPT Production
   - Energy Minimization → Annealing → NVT Equilibration → Production NVT  
   - Repeated **3 times** for statistical reliability
  
  ### Analysis Performed
- From **NPT Production**:
  - Density
  - Box size

- From **NVT Production**:
  - Radial Distribution Function (RDF)
  - Mean Squared Displacement (MSD)
  - Diffusion coefficient (from MSD)
  - Ionic Conductivity
    
## Repository Structure
```
ionic_liquid_Na-Pyr13-TFSI/
│
├── simulation_input/
│   ├── input_file/
|      ├── il.ff
|      ├── c3c1pyrr.zmat
|      ├── Na.zmat
|      ├── ntf2.zmat
|   ├── script/
|      ├── fftool
|      ├── mdp_file/
|         └── ...
|
├── simulation_output/
|   ├── NPT_Production/
|      └── ...   
|   ├── NVT_Production/
|      └── ...

│
├── analysis/
|   ├── script/
|      └── ...
│   ├── RDF/
|      └── ...
|   ├── MSD/
|      └── ...
|   ├── Conductivity/
|      └── ...
│
└── docs/
```
