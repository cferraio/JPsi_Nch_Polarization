#!/bin/sh

homedir=`pwd`
cd ..
cd ..
basedir=`pwd`
cd macros/polFit
#storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir
storagedir=$basedir/Psi/Data

########## INPUTS ##########

for nState in 4 5;do

cp ../../interface/rootIncludes.inc               rootIncludes.inc
cp ../../interface/commonVar_Psi$[nState-3]S.h    commonVar.h
cp ../../interface/ToyMC_Psi$[nState-3]S.h        ToyMC.h

SystID=TotalSyst

JobID=$storagedir/SquaredSystApr8_All

if [ $nState -eq 4 ]
then
ptBinMin=1
ptBinMax=12
fi
if [ $nState -eq 5 ]
then
ptBinMin=1
ptBinMax=5
fi

########################################

cd ${homedir}

touch ChangeTGraph.cc
make

mkdir Systematics
mkdir Systematics/${SystID}

SystDir=Systematics/${SystID}/ChangedTGraph

mkdir ${SystDir}


./ChangeTGraph ${JobID}=JobID1 ${SystID}=SystID ${storagedir}=storagedir ${basedir}=basedir ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${nState}nState
rm ChangeTGraph

done


