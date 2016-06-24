#!/bin/bash
gfortran -o meshgen -O2 src/meshabaqus9.f
./meshgen
cp ./mesh/dbcs ./input/dbcs
cp ./mesh/nodes ./input/nodes
cp ./mesh/elements ./input/elements
cp ./mesh/fbcs ./input/fbcs
cp ./mesh/crackloc ./input/crackloc
cp ./mesh/moduli ./input/moduli
cp ./mesh/DtN_dofs ./input/DtN_dofs
cp ./mesh/gcon ./input/gcon
cp ./mesh/gconu ./input/gconu
cp ./mesh/gconmu ./input/gconmu
cp ./mesh/gconvp ./input/gconvp
echo 'TRANSFER COMPLETE'
