storagedir=/tmp_mnt/scratch/knuenz/Polarization/Upsilon
#macdir=/Users/valentinknuenz/usr/local/workspace/UpsilonPol/macros/polFit
AFSdir=/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit

DataToy=ToyMC

for JobID in Toy0128_AN_Eff105NoDiEffNoRho;do

echo $JobID
cp -r ${storagedir}/${DataToy}/${JobID}/*/Figures/ ${AFSdir}/Figures${DataToy}/${JobID}

#rm -r ${AFSdir}/Figures${DataToy}/${JobID}/Figures/

done
