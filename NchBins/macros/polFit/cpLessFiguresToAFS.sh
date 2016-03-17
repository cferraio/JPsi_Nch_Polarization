storagedir=/scratch/knuenz/Polarization/Upsilon/Data
macdir=/afs/hephy.at/scratch/k/knuenz/CMSSW_4_2_4_patch2/src/UpsilonPol/macros/polFit

DataToy=Data

for JobID in Data_TowardsPRL_Aug11_FinalResults_1Sigma_AlteredPPD_BKGlinPLUSRestSquaredGauss_5nRand_BiasCorrection_1S2S3SAug12_1Sig;do
for nState in 1;do
for rap_ in 1;do
for pT_ in 10;do

echo $JobID
mkdir ${macdir}/Figures${DataToy}/${JobID}
mkdir ${macdir}/Figures${DataToy}/${JobID}/Figures

cp ${storagedir}/${JobID}/Figures/fit_CS_costh_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_CS_phi_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_CS_phith_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_HX_costh_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_HX_phi_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_HX_phith_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_PX_costh_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_PX_phi_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_PX_phith_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/lph_vs_lth_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/lphstar_vs_lthstar_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/ltp_vs_lph_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/ltp_vs_lth_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/lph_vs_lth_${nState}SUps_rap${rap_}_pT${pT_}.jpg ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/lphstar_vs_lthstar_${nState}SUps_rap${rap_}_pT${pT_}.jpg ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/ltp_vs_lph_${nState}SUps_rap${rap_}_pT${pT_}.jpg ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/ltp_vs_lth_${nState}SUps_rap${rap_}_pT${pT_}.jpg ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/lph_vs_lth_${nState}SUps_rap${rap_}_pT${pT_}.C ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/lphstar_vs_lthstar_${nState}SUps_rap${rap_}_pT${pT_}.C ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/ltp_vs_lph_${nState}SUps_rap${rap_}_pT${pT_}.C ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/ltp_vs_lth_${nState}SUps_rap${rap_}_pT${pT_}.C ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/ltilde_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_cosalpha_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_background_rap_test_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_background_mass_test_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_background_pT_test_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_background_phiPX_test_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_background_costhPX_test_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures

cp ${storagedir}/${JobID}/Figures/fit_total_rap_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_total_mass_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_total_pT_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_total_phiPX_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/fit_total_costhPX_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures

cp ${storagedir}/${JobID}/Figures/PosteriorDist_lth_CS_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/PosteriorDist_lph_CS_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/PosteriorDist_ltp_CS_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/PosteriorDist_ltilde_CS_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/PosteriorDist_lth_HX_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/PosteriorDist_lph_HX_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/PosteriorDist_ltp_HX_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/PosteriorDist_ltilde_HX_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/PosteriorDist_lth_PX_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/PosteriorDist_lph_PX_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/PosteriorDist_ltp_PX_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures
cp ${storagedir}/${JobID}/Figures/PosteriorDist_ltilde_PX_${nState}SUps_rap${rap_}_pT${pT_}.pdf ${macdir}/Figures${DataToy}/${JobID}/Figures



done
done
done
done


