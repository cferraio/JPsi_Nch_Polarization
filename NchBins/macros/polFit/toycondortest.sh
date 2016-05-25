genMin=1
genMax=50
ptBinMin=1
ptBinMax=2
cpmBinMin=1
cpmBinMax=10
x=0
JobID=FrameworkIII_19May2016/Sig_frame1scen3_Bkg_frame1scen4

rm blah.log

gen_=${genMin}
while [ $gen_ -le ${genMax} ]
do
pT_=${ptBinMin}
while [ $pT_ -le ${ptBinMax} ]
do
cpm_=${cpmBinMin}
while [ $cpm_ -le ${cpmBinMax} ]
do

FILE=/data/users/ferraioc/Polarization/JPsi/NchBins/ToyMC/$JobID/rap1_pT${pT_}_cpm${cpm_}/Generation${gen_}/results.root 


if [ -f "$FILE" ];
then
   echo ""
else
   echo "File $FILE does not exist" >> toycondortest.log
   echo "$pT_ $cpm_ ${gen_}"
   cp runcondorToys.jdl runtoyscondor.jdl
   redo=$[gen_-1]
   echo "$redo condor-simple.py $pT_ $cpm_" >> runtoyscondor.jdl
   echo "Queue 1" >> runtoyscondor.jdl
#   sleep 1
#   condor_submit runtoyscondor.jdl
   x=$[x+1]
fi

echo $x

cpm_=$[cpm_+1]
done
pT_=$[pT_+1]
done
gen_=$[gen_+1]
done