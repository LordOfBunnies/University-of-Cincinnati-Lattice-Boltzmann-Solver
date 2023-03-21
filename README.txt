The University of Cincinnati Lattice Boltzmann Solver (UCLBS) is an open source entropic lattice Boltzmann method (ELBM) solver based on the works of Karlin et al. from ETH Zurich.  The code was created by Dr. Sean Duncan as part of his PhD dissertation (Open Source, High-Speed Lattice Boltzmann Solver) at the University of Cincinnati (UC).  It uses adaptive mesh refinement (AMR) to adjust the mesh 

The world of lattice Boltzmann work can be somewhat insular, so this was created and made open source as an attempt to demystify some of the difficulties of creating LBM software.

Details will be included below.

The software may or may not be supported long term, but it is intended as a base for others to start from.  No guarantee of functionality is made or represented.

A link to my dissertation is below, detailing the work and systems of UCLBS.

http://rave.ohiolink.edu/etdc/view?acc_num=ucin1668636771464341



OVERVIEW
UCLBS is meant to be an aerospace flow solver.  It is capable of simulating flows up to Mach 2.5, but modification may be needed to startup procedure to achieve initial flows of this velocity.  

ELBM is similar to an LES method.  The entropic estimate, alpha, acts similarly to a sub-grid scale model.  As a result, simulations do not take into account turbulence as it will be generated automatically.  In cases of the boundary layer flow, it is necessary to account for turbulent transition as well.  For example, for Poiseiulle flow, an entrance length for transition to turbulence must be accounted for if studying turbulent pipe flow.



INPUT
UCLBS takes several files as input.  The first is a NAMELIST filed which must be called uclbs.in .  This defines many of the flow parameters, mesh characteristics, and other files to be used.  Within uclbs.in, you can add a number of geometry or CGNS files to the read phase.  The second required file is the Amrex input file, and it is added on the command line after invoking UCLBS.  The Amrex input file defines the grid boundaries, node counts, and several refinement parameters.

All input is metric units.



EXTERNAL CODES
UCLBS uses 4 external codes as part of it, and they need to be installed on the system you're intending to use and pointed to when compiled.  These codes are Amrex 22.06 or higher, CGNS 4.10 or higher, OpenMPI or MPICH 4.0 or higher, and HDF5 1.12 or higher.  In addition, two files in Amrex need to updated or replaced with ones included in this software package, detailed in the Amrex_Replacement_Files folder.  Amrex must have HDF5 enabled, parallel processing, and Fortran enabled.



MESHING
Meshing is automatic in UCLBS, but you can define Boxes of Interest (BOI) where the mesh will be refined.  In addition, you can define a coarse mesh and coarse time to allow the flow time to initialize before increasing the refinement.



GEOMETRY
Geometry in UCLBS is done either through .stl files and through CGNS files.  Any number of geometry files can be used, but the simulations done using .stl files are generally simpler cases with only the edges of the mesh being boundary conditions.  All geometry is considered solid, no-slip walls.  

Be aware of your coarse grid size, the maximum lattice length, and the thickness of geometry pieces.  The node identification system is not perfect when dealing with thin pieces of geometry.

As well, nodes are checked against EVERY triangle when refinement occurs, so increasing the number of triangles can slow down the regrid process.



CGNS
Input can also be done using a CGNS file or a mix of CGNS and .stl files. 

If using a CGNS file for input, several requirements exist.  All surfaces in the volume must be meshed with a triangular mesh.  Boundaries must be labeled, inlet_1 (INLET_1 or Inlet_1 also), outlet_1, freestream_1, etc.  Everything else is taken as solid geometry. Inlet and outlet information is defined in the uclbs.in file.  Additionally, avolume mesh MUST exist, but it does not need to be refined.  The volume mesh is not used, but it generally required for the .cgns file to output properly.

Complex surfaces should be meshed with greater refinement to ensure the details are captured.  The viscous wall systems can account for the geometry, but more detail ensures more information is picked up.



AMREX and AMR
UCLBS uses Amrex to handle all the AMR work.  A system created by Dr. Duncan finds the nodes requiring refinement during regrid steps.  

Currently, the system does not allow regrid on null nodes meaning the refined areas will mostly only get smaller.




COMPILING UCLBS
UCLBS requires gcc 8.* or higher as it needs the std-c17 standards.  Currently, the object files are stored under the Debug folder.  UCLBS was created in Eclipse and the makefiles were generated automatically.  As a result, locations for things like the Amrex and CGNS files will need to be updated.  A shell script has been placed in the top level of this folder to update all the required files, so it will need to be modified for your system.

Eventually, I would like to make this more standardized.




USAGE
UCLBS can be used for many purposes, but it is somewhat limited right now.




TEST CASES
Several test cases are included to demonstrate the possibilities of UCLBS.

1) Pipe Flow
Laminar or turbulent pipe flow case.  Laminar requires a long pipe, and some unstable oscillations are not damped out leading to problems eventually.  For the turbulent case, ensure enough length for transition to turbulence to occur as well as fully developing the flow.

2) Boundary Layer Flow
Simple boundary layer flow.  Other side of the geometry can also be used for a lid driven cavity flow.

3) Supersonic Wedge
10o shock generator at Mach 1.5.  Viscous wedge and inviscid sidewalls.

4) Supersonic Nozzle
Mach 1.5 nozzle, starting from quiescent into 5 m/s air.

