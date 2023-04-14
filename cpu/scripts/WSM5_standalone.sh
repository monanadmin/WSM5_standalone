#!/bin/sh

# Script que compila e executa o codigo WSM5_StandAlone:

# Usage: 
#
#        ./WSM5_standalone.sh  
#

DIRHOME=$PWD/..
SCRIPTS=${DIRHOME}/scripts
DATAOUT=${DIRHOME}/dataout
DATAIN=${DIRHOME}/datain
SRC=${DIRHOME}/src
BIN=${DIRHOME}/bin


cd ${BIN}
rm -f ${BIN}/wsm.x
echo "Compilando"
#make clean
make 
sleep 3


cd ${DATAIN}
echo "Executando"
time ${BIN}/wsm.x

echo "done"
