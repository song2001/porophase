#!/bin/bash
echo "File Name?: (<filename>.inp)"
read TEMPFILE
FILE1=$(echo "${TEMPFILE}.inp")
NODESTART="$(grep -n "*Node" < "./mesh/abaqus/$FILE1" | cut -d : -f 1)"
NODEEND="$(grep -n "*Element," < "./mesh/abaqus/$FILE1" | cut -d : -f 1)"
ELESTART="$(grep -n "*Element," < "./mesh/abaqus/$FILE1" | cut -d : -f 1)"
ELEEND="$(grep -n "*End Part" < "./mesh/abaqus/$FILE1" | cut -d : -f 1)"
NODESTART=$(echo "scale=1; $NODESTART + 1" | bc -l)
NODEEND=$(echo "scale=1; $NODEEND - 1" | bc -l)
ELESTART=$(echo "scale=1; $ELESTART + 1" | bc -l)
ELEEND=$(echo "scale=1; $ELEEND - 1" | bc -l)
NUMNODES=$(echo "scale=1; $NODEEND - $NODESTART + 1" | bc -l)
NUMELES=$(echo "scale=1; $ELEEND - $ELESTART + 1" | bc -l)
echo "Nodes: $NUMNODES"
echo "Elements: $NUMELES"
echo $NUMNODES > ./mesh/nodes_in
echo $NUMELES > ./mesh/elements_in
sed -n "$NODESTART","$NODEEND"p "./mesh/abaqus/$FILE1" >> ./mesh/nodes_in
sed -i -e 's/,/ /g' ./mesh/nodes_in
sed -n "$ELESTART","$ELEEND"p "./mesh/abaqus/$FILE1" >> ./mesh/elements_in
sed -i -e 's/,/ /g' ./mesh/elements_in
