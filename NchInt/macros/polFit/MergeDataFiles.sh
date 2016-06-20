#!/bin/sh

########## INPUTS ##########

nState=4
NSigma=3.00 #needed in 2 decimal accuracy (x.yz)

#JobID=Psi$[nState-3]S_${NSigma}Sigma_11Dec2012
JobID=Psi$[nState-3]S_${NSigma}Sigma_11Dec2012_noRhoFactor

nGenerations=10
MergeFiles=1

if [ $nState -eq 4 ] 
then
rapBinMin=1
rapBinMax=2
ptBinMin=1
ptBinMax=12
fi

if [ $nState -eq 5 ] 
then
rapBinMin=1
rapBinMax=3
ptBinMin=2
ptBinMax=6
fi

nSkipGen=0

########################################

TreeID=Psi$[nState-3]S

##############################################

cd ..
cd ..
basedir=$PWD
cd macros/polFit
#storagedir=`more storagedir`/Data #please define the directory storagedir in the file macros/polFit/storagedir
storagedir=$basedir/Psi/Data

mkdir ${storagedir}
mkdir ${storagedir}/${JobID}

cd ${storagedir}/${JobID}
pwd

rap_=${rapBinMin}
while [ $rap_ -le ${rapBinMax} ]
do
pT_=${ptBinMin}
while [ $pT_ -le ${ptBinMax} ]
do

nGen_=${nSkipGen}
nGen_=$[nGen_+1]
while [ $nGen_ -le ${nGenerations} ]
do


resultfilename=resultsMerged_Psi$[nState-3]S_rap${rap_}_pT${pT_}.root


nFitResultName=results_Fit${nGen_}_Psi$[nState-3]S_rap${rap_}_pT${pT_}.root
echo ${nFitResultName}

if [ $MergeFiles -eq 1 ]
then

if [ $nGen_ -eq 1 ]
then
cp ${nFitResultName} ${resultfilename}
fi

if [ $nGen_ -ge 2 ]
then
mv ${resultfilename} BUFFER_${resultfilename}
hadd -f ${resultfilename} BUFFER_${resultfilename} ${nFitResultName}
rm BUFFER_${resultfilename}
fi

fi


nGen_=$[nGen_+1]
done


FinalDestinationName=results_Psi$[nState-3]S_rap${rap_}_pT${pT_}.root

if [ $MergeFiles -eq 1 ]
then
mv ${FinalDestinationName} tmp_${FinalDestinationName}
mv ${resultfilename} ${FinalDestinationName}
fi

pT_=$[pT_+1]
done
rap_=$[rap_+1]
done


