SIMULATION INSTRUCTION
=============================
A. NPT Production
=============================

A.1. Build Box Simulation
--------

STEP 1. Generates input file for Packmol using fftool
--------

- **files** : [Na.zmat](simulation_input/input_file/Na.zmat) , [c3c1pyrr.zmat](simulation_input/input_file/c3c1pyrr.zmat), [ntf2.zmat](simulation_input/input_file/ntf2.zmat), and [il.ff](simulation_input/input_file/il.ff)
- **script** : [fftool](simulation_input/script/FFTool)
- **command** : replace placeholders with the number of ions*
  ```
  ./fftool <n_Na> Na.zmat <n_pyr13> c3c1pyrr.zmat <n_TFSI> ntf2.zmat -r 1.5
  ```
- **example**
  ```
  ./fftool 19 Na.zmat 171 c3c1pyrr.zmat 190 ntf2.zmat -r 1.5 
  ```
  <details>
  <summary>expected terminal output</summary>
    
  <pre>
    density 1.500 mol/L  volume 420679.7 A^3
    molecule_file      species           nmol force_field      natom nbond source  charge
      Na.zmat          Na+                 19 il.ff                1     0 file   +1.0000
      c3c1pyrr.zmat    c3c1pyrr+          171 il.ff               27    27 file   +1.0000
      ntf2.zmat        tf2N-              190 il.ff               15    14 file   -1.0000
    packmol file
    pack.inp
  </pre>
  
  </details>


- **output file** : `pack.inp` (input for packmol)


STEP 2. Generates initial configurations using Packmol
--------
- **files** : `pack.inp`
- **software** : [Packmol](http://www.ime.unicamp.br/~martinez/packmol/)
- **command** :
  - search module 
    ```
    module avail packmol
    ```

    <details>
    <summary>expected terminal output</summary>
        
      <pre>
    ------------------------------------------------------- /mgpfs/apps/modulefiles -------------------------------------------------------
     nuclear/packmol/20.14.4
  
    If the avail list is too long consider trying:
  
    "module --default avail" or "ml -d av" to just list the default modules.
    "module overview" or "ml ov" to display the number of modules for each name.
  
    Use "module spider" to find all possible modules and extensions.
    Use "module keyword key1 key2 ..." to search for all possible modules matching any of the "keys".
  
      </pre>
      </details>

  - load module
    ```
    module load nuclear/packmol/20.14.4
    ```
  - check module
    ```
    module list
    ```
    <details>
    <summary>expected terminal output</summary>
      
    <pre>
    Currently Loaded Modules:
    1) prun/2.2       3) hwloc/2.7.0   5) libfabric/1.13.0   7) singularity/3.7.1   9) readline/8.2
    2) gnu12/12.2.0   4) ucx/1.14.0    6) openmpi4/4.1.4     8) ohpc               10) nuclear/packmol/20.14.4
    </pre>

    
    </details>
    
    > 💡 Note: On HPC systems, you only need to load the required software (Packmol) **once per session**.  
    > If you log out or start a new session, you need to load it again before running the commands.


  - generates initial configurations
    ```
    packmol < pack.inp
    ```
      <details>
      <summary>expected terminal output</summary>
      
      <pre>
        ################################################################################
        
         PACKMOL - Packing optimization for the automated generation of
         starting configurations for molecular dynamics simulations.
         
                                                                      Version 20.3.5 
        
        ################################################################################
        
          Packmol must be run with: packmol < inputfile.inp 
        
          Userguide at: http://m3g.iqm.unicamp.br/packmol 
        
          Reading input file... (Control-C aborts)
          Seed for random number generator:      1234567
          Output file: simbox.xyz
          Reading coordinate file: ilff/Na_pack.xyz
          Number of independent structures:            1
          The structures are: 
          Structure            1 :Na+ il.ff(           1  atoms)
          Maximum number of GENCAN loops for all molecule packing:          200
          Total number of restrictions:            1
          Distance tolerance:    2.5000000000000000     
          Number of molecules of type            1 :           20
          Total number of atoms:           20
          Total number of molecules:           20
          Number of fixed molecules:            0
          Number of free molecules:           20
          Number of variables:          120
          Total number of fixed atoms:            0
          Maximum internal distance of type            1 :    0.0000000000000000     
          All atoms must be within these coordinates: 
           x: [   -998.39999999999998      ,    1001.6000000000000       ] 
           y: [   -998.39999999999998      ,    1001.6000000000000       ] 
           z: [   -998.39999999999998      ,    1001.6000000000000       ] 
          If the system is larger than this, increase the sidemax parameter. 
        
        ################################################################################
        
          Building initial approximation ... 
        
        ################################################################################
        
          Adjusting initial point to fit the constraints 
        
        --------------------------------------------------------------------------------
        
        --------------------------------------------------------------------------------
        
          Molecules of type:            1
        
          Packing:|0                                                             100%|
                  |*******
        
          Restraint-only function value:    7.6356809804166854E-014
          Maximum violation of the restraints:    1.7918856863805915E-014
        
        --------------------------------------------------------------------------------
        
          Rescaling maximum and minimum coordinates... 
          Computing size of patches... 
          Add fixed molecules to permanent arrays... 
          Reseting center of mass... 
        
        --------------------------------------------------------------------------------
        
          Setting initial trial coordinates ... 
        
        --------------------------------------------------------------------------------
        
        --------------------------------------------------------------------------------
        
          Molecules of type:            1
          Adjusting random positions to fit the constraints. 
          Restraint-only function value:    0.0000000000000000     
          Maximum violation of the restraints:    0.0000000000000000     
        
        ################################################################################
        
          Objective function at initial point:    0.0000000000000000     
        
        ################################################################################
        
          Packing molecules of type:            1
        
        ################################################################################
        
        
          Initial approximation is a solution. Nothing to do. 
        
          Current point written to file: simbox.xyz
        --------------------------------------------------------------------------------
          Packing solved for molecules of type           1
          Objective function value:    0.0000000000000000     
          Maximum violation of target distance:    0.0000000000000000     
          Max. constraint violation:    0.0000000000000000     
        --------------------------------------------------------------------------------
        
        ################################################################################
        
          Packing all molecules together 
        
        ################################################################################
        
        
          Initial approximation is a solution. Nothing to do. 
        
          Solution written to file: simbox.xyz
        
        ################################################################################
        
                                         Success! 
                      Final objective function value: .00000E+00
                      Maximum violation of target distance:   0.000000
                      Maximum violation of the constraints: .00000E+00
        
        --------------------------------------------------------------------------------
        
                      Please cite this work if Packmol was useful: 
        
                   L. Martinez, R. Andrade, E. G. Birgin, J. M. Martinez, 
                 PACKMOL: A package for building initial configurations for
                           molecular dynamics simulations. 
                  Journal of Computational Chemistry, 30:2157-2164,2009.
        
        ################################################################################
        
           Running time:    1.27600040E-03  seconds. 
        
        -------------------------------------------------------------------------------
      </pre>
      </details>
- **output file** : `simbox.xyz`, `c3c1pyrr_pack.xyz`, `ntf2_pack.xyz`, `Na_pack.xyz`

STEP 3. Generates input file for Gromacs
--------
- **files** : `simbox.xyz`, [Na.zmat](simulation_input/input_file/Na.zmat) , [c3c1pyrr.zmat](simulation_input/input_file/c3c1pyrr.zmat), [ntf2.zmat](simulation_input/input_file/ntf2.zmat), and [il.ff](simulation_input/input_file/il.ff)
- **script** : [fftool](simulation_input/script/FFTool)
- **command** : replace placeholders with the number of ions*
  ```
  ./fftool <n_Na> Na.zmat <n_pyr13> c3c1pyrr.zmat <n_TFSI> ntf2.zmat -r 1.5 -g
  ```
- **example**
  ```
  ./fftool 19 Na.zmat 171 c3c1pyrr.zmat 190 ntf2.zmat -r 1.5 -g
  ```
  <details>
    <summary>expected terminal output</summary>
    
  <pre>
  density 1.500 mol/L  volume 420679.7 A^3
  molecule_file      species           nmol force_field      natom nbond source  charge
    Na.zmat          Na+                 19 il.ff                1     0 file   +1.0000
    c3c1pyrr.zmat    c3c1pyrr+          171 il.ff               27    27 file   +1.0000
    ntf2.zmat        tf2N-              190 il.ff               15    14 file   -1.0000
  gromacs files
    run.mdp
    field.top
    config.pdb
  </pre>
  </details>
- **output file** : `run.mdp`, `field.top`, and `config.pdb`

STEP 4. Convert the PDB structure into a GROMACS-compatible GRO file
--------
- **files** : `config.pdb`
- **software** : [GROMACS](https://www.gromacs.org)
- **command** :
  - search module
    ```
    module avail gromacs
    ```
    <details>
    <summary>expected terminal output</summary>
      
    <pre>
    ------------------------------------------------------- /mgpfs/apps/modulefiles -------------------------------------------------------
       bioinformatics/gromacs/2023-plumed-mpi    bioinformatics/gromacs/2023-plumed    bioinformatics/gromacs/2023.3-mpi (D)
    
      Where:
       D:  Default Module
    
    If the avail list is too long consider trying:
    
    "module --default avail" or "ml -d av" to just list the default modules.
    "module overview" or "ml ov" to display the number of modules for each name.
    
    Use "module spider" to find all possible modules and extensions.
    Use "module keyword key1 key2 ..." to search for all possible modules matching any of the "keys".
    </pre>
    </details>
  
  - module load
    ```
    module load bioinformatics/gromacs/2023.3-mpi
    ```
  - module check
    ```
    module list
    ```
    <details>
    <summary>expected terminal output</summary>
      
    <pre>
    Currently Loaded Modules:
      1) prun/2.2       4) ucx/1.14.0         7) singularity/3.7.1  10) mpi/2021.10.0      13) bioinformatics/plumed/2.9.0
      2) gnu12/12.2.0   5) libfabric/1.13.0   8) ohpc               11) fftw/3.3.8-shared  14) bioinformatics/gromacs/2023.3-mpi
      3) hwloc/2.7.0    6) openmpi4/4.1.4     9) python/3.9.16      12) gcc/12.2.0

    </pre>
    </details>
 
    > 💡 Note: On HPC systems, you only need to load the required software (Gromacs) **once per session**.  
    > If you log out or start a new session, you need to load it again before running the commands.

  - Convert the .pdb to .gro format
    ```
    gmx_mpi editconf -f <file.pdb> -o <file.gro>
    ```
- **example**
  ```
  gmx_mpi editconf -f config.pdb -o config.gro
  ```
  <details>
  <summary>expected terminal output</summary>
      
  <pre>
  Command line:
    gmx editconf -f config.pdb -o config.gro
  
  Note that major changes are planned in future for editconf, to improve usability and utility.
  WARNING: all CONECT records are ignored
  Read 7486 atoms
  Volume: 420.678 nm^3, corresponds to roughly 189300 electrons
  No velocities found
  
  GROMACS reminds you: "A C program is like a fast dance on a newly waxed dance floor by people carrying razors." (Waldi Ravens)

  </pre>
  </details>

- **output file** :`config.gro`
---

A.2. Charge Scaling
--------
- **files** : `field.top`
- **script** : [charge_scaling.py](simulation_input/script/charge_scaling)
- **command** : replace placeholders with scaling factor
  ```
  python3 scaling.py field.top <scaling_factor>
  ```
  
- **example**
    ```
    python3 scaling.py field.top 0.70
    ```
   <details>
    <summary>expected terminal output</summary>
      
    <pre>
    Reading input file: field.top
    Scaling charges with factor: 0.7
    Scaling charges completed:
    - Total atoms scaled in [ atoms ]: 43
    - Total atomtypes scaled in [ atomtypes ]: 13
    Writing scaled output to: field_scaled_0.7000.top
    Done!
    </pre>
    
    </details>

  - **output file** : `field_scaled_0.7000.top`
  
  ---

A.3 Energy Minimization
--------
- **software** : [GROMACS](https://www.gromacs.org)
  
STEP 1. Gromacs Preprocessor / Input Preparation
--------
- **files** : `config.gro`, `field_scaled_0.7000.top`, `em.mdp`
- **command** :

  ```
  gmx_mpi grompp -f em.mdp -p field_scaled_0.7000.top -c config.gro -o em.tpr
  ```
- **output file** : `em.tpr`

STEP 2. Running the Energy Minimization
-------
- **files** : `em.tpr`
- **command** :

  ```
  gmx_mpi mdrun -v -deffnm em
  ```
  > 💡 Note: Running GROMACS simulations (using gmx mdrun) can take a long time  
  > For faster and more convenient execution on HPC systems, it is recommended to run it using SLURM.


