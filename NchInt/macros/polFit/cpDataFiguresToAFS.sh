storagedir=/tmp_mnt/scratch/knuenz/Polarization/Upsilon
#macdir=/Users/valentinknuenz/usr/local/workspace/UpsilonPol/macros/polFit
AFSdir=/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit

DataToy=Data

for JobID in Data_TheGreatRun_10B_Apr18_NewestCentrals;do

echo $JobID
cp -r ${storagedir}/${DataToy}/${JobID}/Figures/ ${AFSdir}/Figures${DataToy}/${JobID}

#rm -r ${AFSdir}/Figures${DataToy}/${JobID}/Figures/

done



cp ../Addition1_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap1_pT7.root results_Fit26_1SUps_rap1_pT7.root
cp ../Addition1_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap1_pT7.root results_Fit27_1SUps_rap1_pT7.root
cp ../Addition2_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap1_pT7.root results_Fit28_1SUps_rap1_pT7.root
cp ../Addition2_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap1_pT7.root results_Fit29_1SUps_rap1_pT7.root
cp ../Addition3_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap1_pT7.root results_Fit30_1SUps_rap1_pT7.root
cp ../Addition3_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap1_pT7.root results_Fit31_1SUps_rap1_pT7.root
cp ../Addition4_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap1_pT7.root results_Fit32_1SUps_rap1_pT7.root
cp ../Addition4_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap1_pT7.root results_Fit33_1SUps_rap1_pT7.root
cp ../Addition5_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap1_pT7.root results_Fit34_1SUps_rap1_pT7.root
cp ../Addition5_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap1_pT7.root results_Fit35_1SUps_rap1_pT7.root
cp ../Addition6_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap1_pT7.root results_Fit36_1SUps_rap1_pT7.root
cp ../Addition6_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap1_pT7.root results_Fit37_1SUps_rap1_pT7.root
cp ../Addition7_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap1_pT7.root results_Fit38_1SUps_rap1_pT7.root
cp ../Addition7_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap1_pT7.root results_Fit39_1SUps_rap1_pT7.root
cp ../Addition8_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap1_pT7.root results_Fit40_1SUps_rap1_pT7.root
cp ../Addition8_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap1_pT7.root results_Fit41_1SUps_rap1_pT7.root
cp ../Addition9_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap1_pT7.root results_Fit42_1SUps_rap1_pT7.root
cp ../Addition9_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap1_pT7.root results_Fit43_1SUps_rap1_pT7.root
cp ../Addition10_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap1_pT7.root results_Fit44_1SUps_rap1_pT7.root
cp ../Addition10_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap1_pT7.root results_Fit45_1SUps_rap1_pT7.root
cp ../Addition11_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap1_pT7.root results_Fit46_1SUps_rap1_pT7.root
cp ../Addition11_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap1_pT7.root results_Fit47_1SUps_rap1_pT7.root
cp ../Addition12_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap1_pT7.root results_Fit48_1SUps_rap1_pT7.root
cp ../Addition12_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap1_pT7.root results_Fit49_1SUps_rap1_pT7.root
cp ../Addition13_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap1_pT7.root results_Fit50_1SUps_rap1_pT7.root
cp ../Addition13_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap1_pT7.root results_FitX1_1SUps_rap1_pT7.root
cp ../Addition14_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap1_pT7.root results_FitX2_1SUps_rap1_pT7.root
cp ../Addition14_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap1_pT7.root results_FitX3_1SUps_rap1_pT7.root
cp ../Addition15_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap1_pT7.root results_FitX4_1SUps_rap1_pT7.root
cp ../Addition15_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap1_pT7.root results_FitX5_1SUps_rap1_pT7.root

cp ../Addition4_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap2_pT7.root results_Fit33_1SUps_rap2_pT7.root
cp ../Addition5_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap2_pT7.root results_Fit34_1SUps_rap2_pT7.root
cp ../Addition5_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap2_pT7.root results_Fit35_1SUps_rap2_pT7.root
cp ../Addition6_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap2_pT7.root results_Fit36_1SUps_rap2_pT7.root
cp ../Addition6_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap2_pT7.root results_Fit37_1SUps_rap2_pT7.root
cp ../Addition7_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap2_pT7.root results_Fit38_1SUps_rap2_pT7.root
cp ../Addition7_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap2_pT7.root results_Fit39_1SUps_rap2_pT7.root
cp ../Addition8_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap2_pT7.root results_Fit40_1SUps_rap2_pT7.root
cp ../Addition8_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap2_pT7.root results_Fit41_1SUps_rap2_pT7.root
cp ../Addition9_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap2_pT7.root results_Fit42_1SUps_rap2_pT7.root
cp ../Addition9_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap2_pT7.root results_Fit43_1SUps_rap2_pT7.root
cp ../Addition10_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap2_pT7.root results_Fit44_1SUps_rap2_pT7.root
cp ../Addition10_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap2_pT7.root results_Fit45_1SUps_rap2_pT7.root
cp ../Addition11_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap2_pT7.root results_Fit46_1SUps_rap2_pT7.root
cp ../Addition11_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap2_pT7.root results_Fit47_1SUps_rap2_pT7.root
cp ../Addition12_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap2_pT7.root results_Fit48_1SUps_rap2_pT7.root
cp ../Addition12_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap2_pT7.root results_Fit49_1SUps_rap2_pT7.root
cp ../Addition13_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap2_pT7.root results_Fit50_1SUps_rap2_pT7.root

cp ../Addition13_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap2_pT7.root results_FitX1_1SUps_rap2_pT7.root
cp ../Addition14_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap2_pT7.root results_FitX2_1SUps_rap2_pT7.root
cp ../Addition14_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap2_pT7.root results_FitX3_1SUps_rap2_pT7.root
cp ../Addition1_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap2_pT7.root results_FitX4_1SUps_rap2_pT7.root
cp ../Addition1_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap2_pT7.root results_FitX5_1SUps_rap2_pT7.root
cp ../Addition2_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap2_pT7.root results_FitX6_1SUps_rap2_pT7.root
cp ../Addition2_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap2_pT7.root results_FitX7_1SUps_rap2_pT7.root
cp ../Addition3_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap2_pT7.root results_FitX8_1SUps_rap2_pT7.root
cp ../Addition3_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit2_1SUps_rap2_pT7.root results_FitX9_1SUps_rap2_pT7.root
cp ../Addition4_Ups1S_Data_TheGreatRun_10B_May11_NewestCentrals/results_Fit1_1SUps_rap2_pT7.root results_FitX10_1SUps_rap2_pT7.root



cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit1_1SUps_rap1_pT6.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit37_1SUps_rap1_pT6.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit2_1SUps_rap1_pT6.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit38_1SUps_rap1_pT6.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit3_1SUps_rap1_pT6.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit39_1SUps_rap1_pT6.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit4_1SUps_rap1_pT6.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit40_1SUps_rap1_pT6.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit5_1SUps_rap1_pT6.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit41_1SUps_rap1_pT6.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit6_1SUps_rap1_pT6.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit42_1SUps_rap1_pT6.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit7_1SUps_rap1_pT6.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit43_1SUps_rap1_pT6.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit8_1SUps_rap1_pT6.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit44_1SUps_rap1_pT6.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit9_1SUps_rap1_pT6.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit45_1SUps_rap1_pT6.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit10_1SUps_rap1_pT6.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit46_1SUps_rap1_pT6.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit11_1SUps_rap1_pT6.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit47_1SUps_rap1_pT6.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit12_1SUps_rap1_pT6.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit48_1SUps_rap1_pT6.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit13_1SUps_rap1_pT6.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit49_1SUps_rap1_pT6.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit14_1SUps_rap1_pT6.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit50_1SUps_rap1_pT6.root




cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit1_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit37_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit2_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit38_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit3_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit39_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit4_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit40_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit5_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit41_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit6_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit42_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit7_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit43_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit8_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit44_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit9_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit45_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit10_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit46_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit11_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit47_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit12_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit48_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit13_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit49_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit14_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit50_1SUps_rap1_pT7.root

cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit15_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit29_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit16_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit30_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit17_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit31_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit18_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit32_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit19_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit33_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit20_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit34_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit21_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit35_1SUps_rap1_pT7.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit22_1SUps_rap1_pT7.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit36_1SUps_rap1_pT7.root

cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit1_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit29_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit2_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit30_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit3_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit31_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit4_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit32_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit5_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit33_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit6_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit34_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit7_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit35_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit8_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit36_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit9_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit37_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit10_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit38_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit11_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit39_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit12_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit40_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit13_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit41_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit14_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit42_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit15_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit43_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit16_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit44_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit17_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit45_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit18_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit46_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit19_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit47_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit20_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit48_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit21_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit49_2SUps_rap1_pT10.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit22_2SUps_rap1_pT10.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit50_2SUps_rap1_pT10.root




cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit1_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit27_2SUps_rap1_pT9.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit2_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit28_2SUps_rap1_pT9.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit3_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit29_2SUps_rap1_pT9.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit4_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit30_2SUps_rap1_pT9.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit5_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit31_2SUps_rap1_pT9.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit6_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit32_2SUps_rap1_pT9.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit7_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit33_2SUps_rap1_pT9.root
cp Additions_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit8_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit34_2SUps_rap1_pT9.root
cp Additions1_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit1_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit35_2SUps_rap1_pT9.root
cp Additions1_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit2_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit36_2SUps_rap1_pT9.root
cp Additions1_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit3_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit37_2SUps_rap1_pT9.root
cp Additions1_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit4_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit38_2SUps_rap1_pT9.root
cp Additions1_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit5_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit39_2SUps_rap1_pT9.root
cp Additions1_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit6_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit40_2SUps_rap1_pT9.root
cp Additions1_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit7_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit41_2SUps_rap1_pT9.root
cp Additions1_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit8_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit42_2SUps_rap1_pT9.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit1_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit43_2SUps_rap1_pT9.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit2_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit44_2SUps_rap1_pT9.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit3_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit45_2SUps_rap1_pT9.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit4_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit46_2SUps_rap1_pT9.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit5_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit47_2SUps_rap1_pT9.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit6_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit48_2SUps_rap1_pT9.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit7_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit49_2SUps_rap1_pT9.root
cp Additions2_Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit8_2SUps_rap1_pT9.root Data_TheGreatApproval_June17_NewestCentrals_SOFT_SOFT/results_Fit50_2SUps_rap1_pT9.root

scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_1SUps_rap1_pT6.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_1SUps_rap1_pT7.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_1SUps_rap1_pT8.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_1SUps_rap1_pT9.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_1SUps_rap1_pT10.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_1SUps_rap2_pT6.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_1SUps_rap2_pT7.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_1SUps_rap2_pT8.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_1SUps_rap2_pT9.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_1SUps_rap2_pT10.root .

scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_2SUps_rap1_pT6.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_2SUps_rap1_pT7.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_2SUps_rap1_pT8.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_2SUps_rap1_pT9.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_2SUps_rap1_pT10.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_2SUps_rap2_pT6.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_2SUps_rap2_pT7.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_2SUps_rap2_pT8.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_2SUps_rap2_pT9.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_2SUps_rap2_pT10.root .

scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_3SUps_rap1_pT6.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_3SUps_rap1_pT7.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_3SUps_rap1_pT8.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_3SUps_rap1_pT9.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_3SUps_rap1_pT10.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_3SUps_rap2_pT6.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_3SUps_rap2_pT7.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_3SUps_rap2_pT8.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_3SUps_rap2_pT9.root .
scp knuenz@lxplus250.cern.ch:/tmp/zhlinl/rootfiles/results_3SUps_rap2_pT10.root .


rfcp /castor/cern.ch/user/h/hwoehri/Polarization/UpsGun1S_RealisticPT/MergedTTrees_Aug2012/onia2MuMu_tree_Ups1S_9p5-12p5_merged_4.root .
mv onia2MuMu_tree_Ups1S_* /scratch/knuenz/Polarization/RootInput/Upsilon/
cd /scratch/knuenz/Polarization/RootInput/Upsilon/
hadd -f onia2MuMu_tree_Ups1S_pT6.root onia2MuMu_tree_Ups1S_9p5-12p5_merged_2.root onia2MuMu_tree_Ups1S_9p5-12p5_merged_3.root onia2MuMu_tree_Ups1S_9p5-12p5_merged_5.root onia2MuMu_tree_Ups1S_9p5-12p5_merged_6.root onia2MuMu_tree_Ups1S_9p5-12p5_merged_7.root 

prep6

cp DataFiles/SetOfCuts11_MCclosure_Ups1S_Aug12_pT6/tmpFiles/data_Ups_rap1_pT6.root DataFiles/SetOfCuts11_MCclosure_Ups1S_Aug12/tmpFiles/data_Ups_rap1_pT6.root
cp DataFiles/SetOfCuts11_MCclosure_Ups1S_Aug12_pT6/tmpFiles/data_Ups_rap2_pT6.root DataFiles/SetOfCuts11_MCclosure_Ups1S_Aug12/tmpFiles/data_Ups_rap2_pT6.root

Trim 1 3 10
check
run



sh PlotDataFits.sh
sh EvalSyst.sh
mv MCclosure_Aug11_Ups3S_10MCSig_TO_MCclosure_Aug11_Ups3S_1DataSig/ BiasCorrectionAug12_3S_1Sig
mv MCclosure_Aug11_Ups3S_10MCSig_TO_MCclosure_Aug11_Ups3S_3DataSig/ BiasCorrectionAug12_3S_3Sig
cp BiasCorrectionAug12_3S_1Sig/TGraphResults_3SUps.root BiasCorrectionAug12_3S_1Sig/TGraphResults_1SUps.root
cp BiasCorrectionAug12_3S_1Sig/TGraphResults_3SUps.root BiasCorrectionAug12_3S_1Sig/TGraphResults_2SUps.root
cp BiasCorrectionAug12_3S_3Sig/TGraphResults_3SUps.root BiasCorrectionAug12_3S_3Sig/TGraphResults_1SUps.root
cp BiasCorrectionAug12_3S_3Sig/TGraphResults_3SUps.root BiasCorrectionAug12_3S_3Sig/TGraphResults_2SUps.root

sh AverageSyst.sh
mv AverageSyst BiasCorrection_1S2S3SAug12_1Sig
mv AverageSyst BiasCorrection_1S2S3SAug12_3Sig
cp -r BiasCorrection_1S2S3SAug12_1Sig /scratch/knuenz/Polarization/Upsilon/Data/
cp -r BiasCorrection_1S2S3SAug12_3Sig /scratch/knuenz/Polarization/Upsilon/Data/

Correct results EvalSyst.sh (after change)

cp -r Systematics/BiasCorrectedResults/Data_TowardsPRL_Aug11_FinalResults_1Sigma_TO_BiasCorrection_1S2S3SAug12_1Sig/ /scratch/knuenz/Polarization/Upsilon/Data/
cp -r Systematics/BiasCorrectedResults/Data_TowardsPRL_Aug11_FinalResults_3Sigma_TO_BiasCorrection_1S2S3SAug12_3Sig/ /scratch/knuenz/Polarization/Upsilon/Data/

sh PlotResults.sh





mv 'Photo on 2011-08-22 at 17.54.jpg' photo1.jpg
mv 'Photo on 2012-01-07 at 20.12.jpg' photo2.jpg
mv 'Foto am 04-07-2010 um 19.12 #2.jpg' photo3.jpg
mv 'Foto am 04-07-2010 um 19.12 #3.jpg' photo4.jpg
mv 'Foto am 04-07-2010 um 19.12.jpg' photo5.jpg
mv 'Foto am 04-07-2010 um 19.13.jpg' photo6.jpg
mv 'Foto am 04-07-2010 um 19.26.jpg' photo7.jpg
mv 'Foto am 11-05-2010 um 23.46.jpg' photo8.jpg
mv 'Foto am 11-05-2010 um 23.48.jpg' photo9.jpg
mv 'Foto am 26-05-2010 um 22.39 #2.jpg' photo10.jpg
mv 'Foto am 26-05-2010 um 22.54 #2.jpg' photo11.jpg
mv 'Foto am 26-05-2010 um 22.55.jpg' photo12.jpg
mv 'Foto am 26-05-2010 um 22.57 #2.jpg' photo13.jpg
mv 'Foto am 26-05-2010 um 22.57 #4.jpg' photo14.jpg
mv 'Foto am 26-05-2010 um 22.58 #2.jpg' photo15.jpg
mv 'Foto am 26-05-2010 um 22.58 #3.jpg' photo16.jpg
mv 'Foto am 26-05-2010 um 22.59 #2.jpg' photo17.jpg
mv 'Foto am 26-05-2010 um 22.59.jpg' photo18.jpg
mv '04-button1.JPG' photo19.jpg
mv '05-button2.JPG' photo20.jpg
mv '06-button3.JPG' photo21.jpg
mv '06-button4.JPG' photo22.jpg
mv '14-bokuball.jpg' photo23.jpg
mv '4 am 11-05-2010 um 23.46 #2.jpg' photo24.jpg
mv '4 am 11-05-2010 um 23.46 #3.jpg' photo25.jpg
mv '4 am 11-05-2010 um 23.46 #4.jpg' photo26.jpg
mv '4 am 11-05-2010 um 23.46.jpg' photo27.jpg
mv '4 am 15-04-2011 um 21.54 #2.jpg' photo28.jpg
mv '4 am 15-04-2011 um 21.54 #3.jpg' photo29.jpg
mv '4 am 15-04-2011 um 21.54 #4.jpg' photo30.jpg
mv '4 am 15-04-2011 um 21.54 #5.jpg' photo31.jpg
mv '4 am 15-04-2011 um 21.54 #6.jpg' photo32.jpg
mv '4 am 15-04-2011 um 21.54 #7.jpg' photo33.jpg
mv '4 am 15-04-2011 um 21.54 #8.jpg' photo34.jpg
mv '4 am 15-04-2011 um 21.54.jpg' photo35.jpg
mv '4 am 15-04-2011 um 21.55 #10.jpg' photo36.jpg
mv '4 am 15-04-2011 um 21.55 #11.jpg' photo37.jpg
mv '4 am 15-04-2011 um 21.55 #12.jpg' photo38.jpg
mv '4 am 15-04-2011 um 21.55 #2.jpg' photo39.jpg
mv '4 am 15-04-2011 um 21.55 #3.jpg' photo40.jpg
mv '4 am 15-04-2011 um 21.55 #4.jpg' photo41.jpg
mv '4 am 15-04-2011 um 21.55 #5.jpg' photo42.jpg
mv '4 am 15-04-2011 um 21.55 #6.jpg' photo43.jpg
mv '4 am 15-04-2011 um 21.55 #7.jpg' photo44.jpg
mv '4 am 15-04-2011 um 21.55 #8.jpg' photo45.jpg
mv '4 am 15-04-2011 um 21.55 #9.jpg' photo46.jpg
mv '4 am 15-04-2011 um 21.55.jpg' photo47.jpg
mv '4 am 23-07-2010 um 18.19 #2.jpg' photo48.jpg
mv '4 am 23-07-2010 um 18.19.jpg' photo49.jpg
mv '4 am 29-01-2011 um 13.49 #2.jpg' photo50.jpg
mv '4 am 29-01-2011 um 13.49 #3.jpg' photo51.jpg
mv '4 am 29-01-2011 um 13.49 #4.jpg' photo52.jpg
mv '4 am 29-01-2011 um 13.49.jpg' photo53.jpg
mv '4 am 29-01-2011 um 13.50 #10.jpg' photo54.jpg
mv '4 am 29-01-2011 um 13.50 #11.jpg' photo55.jpg
mv '4 am 29-01-2011 um 13.50 #12.jpg' photo56.jpg
mv '4 am 29-01-2011 um 13.50 #2.jpg' photo57.jpg
mv '4 am 29-01-2011 um 13.50 #3.jpg' photo58.jpg
mv '4 am 29-01-2011 um 13.50 #4.jpg' photo59.jpg
mv '4 am 29-01-2011 um 13.50 #5.jpg' photo60.jpg
mv '4 am 29-01-2011 um 13.50 #6.jpg' photo61.jpg
mv '4 am 29-01-2011 um 13.50 #7.jpg' photo62.jpg
mv '4 am 29-01-2011 um 13.50 #8.jpg' photo63.jpg
mv '4 am 29-01-2011 um 13.50 #9.jpg' photo64.jpg
mv '4 am 29-01-2011 um 13.50.jpg' photo65.jpg
mv '4-up on 2012-04-15 at 21.57 #2.jpg' photo66.jpg
mv '4-up on 2012-04-15 at 21.57 #3.jpg' photo67.jpg
mv '4-up on 2012-04-15 at 21.57 #4.jpg' photo68.jpg
mv '4-up on 2012-04-15 at 21.57 #5.jpg' photo69.jpg
mv '4-up on 2012-04-15 at 21.57 #6.jpg' photo70.jpg
mv '4-up on 2012-04-15 at 21.57 #7.jpg' photo71.jpg
mv '4-up on 2012-04-15 at 21.57 #8.jpg' photo72.jpg
mv '4-up on 2012-04-15 at 21.57.jpg' photo73.jpg
mv '4-up on 2012-04-15 at 21.58 #2.jpg' photo74.jpg
mv '4-up on 2012-04-15 at 21.58 #3.jpg' photo75.jpg
mv '4-up on 2012-04-15 at 21.58 #4.jpg' photo76.jpg
mv '4-up on 2012-04-15 at 21.58 #5.jpg' photo77.jpg
mv '4-up on 2012-04-15 at 21.58 #6.jpg' photo78.jpg
mv '4-up on 2012-04-15 at 21.58 #7.jpg' photo79.jpg
mv '4-up on 2012-04-15 at 21.58 #8.jpg' photo80.jpg
mv '4-up on 2012-04-15 at 21.58.jpg' photo81.jpg
mv 'Bildschirmfoto 2010-08-10 um 19.47.png' photo82.jpg
mv 'Bildschirmfoto 2010-08-11 um 17.39.16.png' photo83.jpg
mv 'Bildschirmfoto 2010-08-15 um 18.27.22.png' photo84.jpg
mv 'Bildschirmfoto 2010-08-16 um 18.40.58.png' photo85.jpg
mv 'Bildschirmfoto 2010-08-18 um 16.15.32 2.png' photo86.jpg
                                                      mv 'IMG_0229 2.JPG' IMG_0229_2.JPG
mv 'IMG_0232 2.JPG' IMG_0232_2.JPG
mv 'IMG_0233 2.JPG' IMG_0233_2.JPG
mv 'IMG_0234 2.JPG' IMG_0234_2.JPG
mv 'IMG_0235 2.JPG' IMG_0235_2.JPG
mv 'IMG_0277 2.JPG' IMG_0277_2.JPG
mv 'IMG_0282 2.JPG' IMG_0282_2.JPG
mv 'IMG_0302 2.JPG' IMG_0302_2.JPG
mv 'IMG_0310 2.JPG' IMG_0310_2.JPG
mv 'IMG_0311 2.JPG' IMG_0311_2.JPG
mv 'IMG_0380 2.JPG' IMG_0380_2.JPG
mv 'IMG_0412 2.JPG' IMG_0412_2.JPG
mv 'IMG_0413 2.JPG' IMG_0413_2.JPG
mv 'IMG_0414 2.JPG' IMG_0414_2.JPG
mv 'IMG_0415 2.JPG' IMG_0415_2.JPG
mv 'IMG_0416 2.JPG' IMG_0416_2.JPG
mv 'IMG_0417 2.JPG' IMG_0417_2.JPG
mv 'IMG_0418 2.JPG' IMG_0418_2.JPG
mv 'IMG_0419 2.JPG' IMG_0419_2.JPG
mv 'IMG_0505 2.JPG' IMG_0505_2.JPG
mv 'IMG_0530 2.JPG' IMG_0530_2.JPG
mv 'IMG_0794 2.JPG' IMG_0794_2.JPG
mv 'IMG_0797 2.JPG' IMG_0797_2.JPG
mv 'IMG_0803 2.JPG' IMG_0803_2.JPG
mv 'IMG_0844 2.JPG' IMG_0844_2.JPG
mv 'IMG_0882 2.JPG' IMG_0882_2.JPG
mv 'IMG_0902 2.JPG' IMG_0902_2.JPG
mv 'IMG_0917 2.JPG' IMG_0917_2.JPG
mv 'IMG_1049 2.JPG' IMG_1049_2.JPG
mv 'IMG_1337 2.JPG' IMG_1337_2.JPG
mv 'IMG_1371 2.JPG' IMG_1371_2.JPG
mv 'IMG_1381 2.JPG' IMG_1381_2.JPG
mv 'IMG_1496 2.JPG' IMG_1496_2.JPG
mv 'IMG_1503 2.JPG' IMG_1503_2.JPG                                                     