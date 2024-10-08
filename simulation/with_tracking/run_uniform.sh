#!/usr/bin/bash

#N.B. STAR version should be DEV

nevents=${1:-10}
echo "Running $nevents events per setting."

#Set particle and associated running tag
particle=${2:-mu_minus}
echo "Running single ${particle} simulation."

echo "Running with E = [1,100] GeV and Eta = [+2.4,+4.5]!"

if [ "${particle}" = "mu_minus" ]; then
    ptag=mu-
elif [ "${particle}" = "pi_minus" ]; then
    ptag=pi-
fi

echo "Particle tag is ${ptag}."
echo ""
sleep 5

#Single particle simulation ( vertex = (0,0,0) )
run_geant=1
if [ "$run_geant" -eq 1 ] ; then
	echo ""
	echo "Starting Geant simulation"
	root4star -b -q runSimUniform.C\(${nevents},1,\"${ptag}\",1\)
else
	echo ""
	echo "Linking to existing Geant simulation file"
	ln -s output/${particle}/uniform/${ptag}.run1.fzd
fi

#Ideal ftt seeding
run_ideal_ftt=0
if [ "$run_ideal_ftt" -eq 1 ] ; then
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
	mv StFwdFitQAMaker.root ${output_dir} #Created by StFwdFitQAMaker
	mv fcstrk.root ${output_dir}  #Created by StFcsTrackMatchMaker
	mv fwdHists.root ${output_dir} #Created if TrackMaker creates histograms
	mv fwdtree.root ${output_dir} #Created if TrackMaker creates TTree
	mv mu-.hist.root ${output_dir}
        mv mu-.run1.geom.root ${output_dir}
        mv mu-.run1.picoDst.root ${output_dir} #Created by PicoDSTMaker
        mv mu-.runco.root ${output_dir}
fi

#Ideal fst seeding
run_ideal_fst=0
if [ "$run_ideal_fst" -eq 1 ] ; then
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
	mv StFwdFitQAMaker.root ${output_dir} #Created by StFwdFitQAMaker
	mv fcstrk.root ${output_dir}  #Created by StFcsTrackMatchMaker
	mv fwdHists.root ${output_dir} #Created if TrackMaker creates histograms
	mv fwdtree.root ${output_dir} #Created if TrackMaker creates TTree
	mv mu-.hist.root ${output_dir}
        mv mu-.run1.geom.root ${output_dir}
        mv mu-.run1.picoDst.root ${output_dir} #Created by PicoDSTMaker
        mv mu-.runco.root ${output_dir}
fi

#Realistic ftt seeding
run_real_ftt=0
if [ "$run_real_ftt" -eq 1 ] ; then
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
	mv StFwdFitQAMaker.root ${output_dir} #Created by StFwdFitQAMaker
	mv fcstrk.root ${output_dir}  #Created by StFcsTrackMatchMaker
	mv fwdHists.root ${output_dir} #Created if TrackMaker creates histograms
	mv fwdtree.root ${output_dir} #Created if TrackMaker creates TTree
	mv mu-.hist.root ${output_dir}
        mv mu-.run1.geom.root ${output_dir}
        mv mu-.run1.picoDst.root ${output_dir} #Created by PicoDSTMaker
        mv mu-.runco.root ${output_dir}
fi

#Realistic fst seeding
run_real_fst=1
if [ "$run_real_fst" -eq 1 ] ; then
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
	mv StFwdFitQAMaker.root ${output_dir} #Created by StFwdFitQAMaker
	mv fcstrk.root ${output_dir}  #Created by StFcsTrackMatchMaker
	mv fwdHists.root ${output_dir} #Created if TrackMaker creates histograms
	mv fwdtree.root ${output_dir} #Created if TrackMaker creates TTree
	mv mu-.hist.root ${output_dir}
	mv mu-.run1.geom.root ${output_dir}
	mv mu-.run1.picoDst.root ${output_dir} #Created by PicoDSTMaker
	mv mu-.runco.root ${output_dir}
fi

#Move generated and Geant files
if [ "$run_geant" -eq 1 ] ; then
	mv ${ptag}.run1.fzd output/${particle}/uniform/
	mv ${ptag}.run1.root output/${particle}/uniform/
else
	unlink ${ptag}.run1.fzd
fi

#Move memory log
#mv mem.dat log/

#Delete remaining ROOT files
#rm -f ${ptag}.*.root

echo ""
echo "Done!"

