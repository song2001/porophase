# porophase
FEM code for phase-field fracture of poroelastic media

Must have Petsc installed to run the code
 PETSC_DIR=(local petsc directory)

To start:

bash Gen_Directories.sh to generate output directories
bash Abaqus2fortran.sh and type in one of the Abaqus .inp files
bash MoveInputs.sh to generate input files

complie code by typing:
make

To run code:
bash RunPoro.sh
