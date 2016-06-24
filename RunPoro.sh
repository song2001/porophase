#!/bin/bash
read -r -p "Store Copies? [y/N] " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
then
  cp -r ./paraview1 ./temp
  cp -r ./results1 ./temp
  cp -r ./paraview2 ./temp
  cp -r ./results2 ./temp
  cp -r ./paraview3 ./temp
  cp -r ./results3 ./temp
else
  echo 'No Copies Stored'
fi
read -r -p "Clear Output? [y/N] " response
if [[ $response =~ ^([yY][eE][sS]|[yY])$ ]]
then
  echo 'CLEARING OUTPUT...'
  rm paraview1/*
  rm results1/sol/*
  rm paraview2/*
  rm results2/sol/*
  rm paraview3/*
  rm results3/sol/*
else
  echo 'Output Not Cleared'
fi
#$FC -O2 -o poro src/porostokes.f
#echo
#echo 'COMPILATION COMPLETE'
#echo
read -r -p "Press any key to run (press 'q' to quit) " response
if [[ $response != 'q' ]]
then
  mpiexec -n 1 ./testp -pc_type lu -ksp_type preonly -pc_factor_mat_solver_package mumps
#./poro
fi
