#!/usr/bin/bash

#N.B. STAR version should be DEV

nevents=${1:-10}
echo "Running $nevents events per setting."

#Set particle and associated running tag
particle=${2:-mu_minus}
echo "Running single ${particle} simulation."

echo "Running with E = [1,100] GeV and Eta = [+2.25,+4.75]!"

if [ "${particle}" = "mu_minus" ]; then
    ptag=mu-
fi

echo "Particle tag is ${ptag}."
echo ""

#Single particle simulation ( vertex = (0,0,0) )
echo ""
echo "Starting Geant simulation"
root4star -b -q runSimUniform.C\(${nevents},1,\"${ptag}\",1\)

#Ideal ftt seeding
echo ""
echo "Starting ideal FTT seeded reconstruction"
root4star -b -q 'recon_uniform.C('${nevents}',"StFwdTrackMaker_ideal_sim_ftt_seed.root", false, true, false,1,"'${ptag}'")' \
| tee log/output_${ptag}_uniform_fttideal.dat

#Create simple ROOT file for analysis
echo ""
echo "Creating ROOT file for analysis"
root4star -b -q 'readMudst.C(0,1,"'${ptag}'.MuDst.root")'

#Move reconstruction files
output_dir=output/${particle}/uniform/ideal_ftt_seed/
mv ${ptag}.MuDst.root ${output_dir} #MuDST file
mv SimpleTree_mudst.root ${output_dir} #Created by my Maker
mv ${ptag}.event.root ${output_dir} #Event file
mv StFwdAnalysisMaker.root ${output_dir} #Created by StFwdAnalysisMaker
mv StFwdFitQAMaker.root ${output_dir} #Created by StFwdFitQAMaker
#mv fcstrk.root ${output_dir}  #Created by StFcsTrackMatchMaker
#mv StFwdTrackMaker_ideal_sim_ftt_seed.root ${output_dir} #Created if TrackMaker creates histograms
#mv fwdtree.root ${output_dir} #Created if TrackMaker creates TTree

#Ideal fst seeding
echo ""
echo "Starting ideal FST seeded reconstruction"
root4star -b -q 'recon_uniform.C('${nevents}',"StFwdTrackMaker_ideal_sim_fst_seed.root", true, true, false,1,"'${ptag}'")' \
| tee log/output_${ptag}_uniform_fstideal.dat

#Create simple ROOT file for analysis
echo ""
echo "Creating ROOT file for analysis"
root4star -b -q 'readMudst.C(0,1,"'${ptag}'.MuDst.root")'

#Move reconstruction files
output_dir=output/${particle}/uniform/ideal_fst_seed/
mv ${ptag}.MuDst.root ${output_dir} #MuDST file
mv SimpleTree_mudst.root ${output_dir} #Created by my Maker
mv ${ptag}.event.root ${output_dir} #Event file
mv StFwdAnalysisMaker.root ${output_dir} #Created by StFwdAnalysisMaker
mv StFwdFitQAMaker.root ${output_dir} #Created by StFwdFitQAMaker
#mv fcstrk.root ${output_dir}  #Created by StFcsTrackMatchMaker
#mv StFwdTrackMaker_ideal_sim_fst_seed.root ${output_dir} #Created if TrackMaker creates histograms
#mv fwdtree.root ${output_dir} #Created if TrackMaker creates TTree

#Realistic ftt seeding
echo ""
echo "Starting realistic FTT seeded reconstruction"
root4star -b -q 'recon_uniform.C('${nevents}',"StFwdTrackMaker_real_sim_ftt_seed.root", false, true, true,1,"'${ptag}'")' \
| tee log/output_${ptag}_uniform_fttreal.dat

#Create simple ROOT file for analysis
echo ""
echo "Creating ROOT file for analysis"
root4star -b -q 'readMudst.C(0,1,"'${ptag}'.MuDst.root")'

#Move reconstruction files
output_dir=output/${particle}/uniform/real_ftt_seed/
mv ${ptag}.MuDst.root ${output_dir} #MuDST file
mv SimpleTree_mudst.root ${output_dir} #Created by my Maker
mv ${ptag}.event.root ${output_dir} #Event file
mv StFwdAnalysisMaker.root ${output_dir} #Created by StFwdAnalysisMaker
mv StFwdFitQAMaker.root ${output_dir} #Created by StFwdFitQAMaker
#mv fcstrk.root ${output_dir}  #Created by StFcsTrackMatchMaker
#mv StFwdTrackMaker_real_sim_ftt_seed.root ${output_dir} #Created if TrackMaker creates histograms
#mv fwdtree.root ${output_dir} #Created if TrackMaker creates TTree

#Realistic fst seeding
echo ""
echo "Starting realistic FST seeded reconstruction"
root4star -b -q 'recon_uniform.C('${nevents}',"StFwdTrackMaker_real_sim_ftt_seed.root", true, true, true,1,"'${ptag}'")' \
| tee log/output_${ptag}_uniform_fstreal.dat

#Create simple ROOT file for analysis
echo ""
echo "Creating ROOT file for analysis"
root4star -b -q 'readMudst.C(0,1,"'${ptag}'.MuDst.root")'

#Move reconstruction files
output_dir=output/${particle}/uniform/real_fst_seed/
mv ${ptag}.MuDst.root ${output_dir} #MuDST file
mv SimpleTree_mudst.root ${output_dir} #Created by my Maker
mv ${ptag}.event.root ${output_dir} #Event file
mv StFwdAnalysisMaker.root ${output_dir} #Created by StFwdAnalysisMaker
mv StFwdFitQAMaker.root ${output_dir} #Created by StFwdFitQAMaker
#mv fcstrk.root ${output_dir}  #Created by StFcsTrackMatchMaker
#mv StFwdTrackMaker_real_sim_fst_seed.root ${output_dir} #Created if TrackMaker creates histograms
#mv fwdtree.root ${output_dir} #Created if TrackMaker creates TTree

#Move generated and Geant files
mv ${ptag}.run1.fzd output/${particle}/uniform/
mv ${ptag}.run1.root output/${particle}/uniform/

#Move memory log
mv mem.dat log/

#Delete remaining ROOT files
rm -f ${ptag}.*.root

echo ""
echo "Done!"

