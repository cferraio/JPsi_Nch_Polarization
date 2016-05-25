fitMin=1
fitMax=25
ptBinMin=1
ptBinMax=2
cpmBinMin=1
cpmBinMax=10
x=0
JobID=19May16_MassUpdateFixedErrBars_FracL25

rm blah.log

fit_=${fitMin}
while [ $fit_ -le ${fitMax} ]
do
pT_=${ptBinMin}
while [ $pT_ -le ${ptBinMax} ]
do
cpm_=${cpmBinMin}
while [ $cpm_ -le ${cpmBinMax} ]
do

FILE=/data/users/ferraioc/Polarization/JPsi/NchBins/$JobID/results_Fit${fit_}_Psi1S_rap1_pT${pT_}_cpm${cpm_}.root

if [ -f "$FILE" ];
then
   echo ""
else
   echo "File $FILE does not exist" >> condortest.log
   echo "${fit_} condor-simple.py $pT_ $cpm_"
   cp runcondorFits.jdl runfitscondor.jdl
   redo=$[fit_-1]
   echo "$redo condor-simple.py $pT_ $cpm_" >> runfitscondor.jdl
   echo "Queue 1" >> runfitscondor.jdl
   sleep 1
   condor_submit runfitscondor.jdl
   x=$[x+1]
fi

echo $x

cpm_=$[cpm_+1]
done
pT_=$[pT_+1]
done
fit_=$[fit_+1]
done