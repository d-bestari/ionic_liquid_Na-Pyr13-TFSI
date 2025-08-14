# MOLECULAR DYNAMICS SIMULATION OF THE EFFECT OF NaTFSI CONCENTRATION AND TEMPERATURE ON IONIC LIQUID Pyr13TFSI 
This repository contains molecular dynamics (MD) simulation data and analysis 
for ionic liquids Na-TFSI and Pyr13-TFSI with varying salt concentrations and temperatures.
The simulations were performed using **GROMACS** with CL&Pol force field parameters.

## Computational Resources
Simulations were performed on the High Performance Computing (HPC) at BRIN (Badan Riset dan Inovasi Nasional), Indonesia.

## Requirements


## Simulation Overview
- **Ionic liquids**: Na-TFSI and Pyr13-TFSI  
- **Salt concentration variations**: 0.0, 0.1, 0.2, 0.3, 0.4 mole fraction  
- **Temperatures**: 300 K, 350 K, 400 K

### Simulation Workflow
1. **Box Volume Determination**
   - Energy minimization
   - Annealing
   - NVT Equilibration
   - NPT Equilibration
   - NPT Equilibration cont.
   - NPT Production
2. **Production NVT**
   - Uses volume obtained from Production NPT
   - EM → Annealing → Equilibration NVT → Production NVT  
   - Repeated **3 times** for statistical reliability
  
  ### Analysis Performed
- From **NPT Production**:
  - Density calculation

- From **NVT Production**:
  - Radial Distribution Function (RDF)
  - Mean Squared Displacement (MSD)
  - Diffusion coefficient (from MSD)
  - Ionic conductivity
