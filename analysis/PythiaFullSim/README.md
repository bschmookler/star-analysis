# Analysis of Pythia Full Simulation results
The scripts in this directory analyze simulated Pythia events which have been passed through the STAR Forward Geant3 model and reconstruction chain. We analyze the MuDst files produced [here](https://github.com/jdbrice/star-sw-1/wiki/Pythia-&-Simulation-samples-and-Productions), which include both Forward track and calorimeter reconstruction.

Analysis codes in this directory
--------------------------------
All analysis codes in this directory use as input the MuDst files converted into a simple TTree. The conversion of the MuDst files into a simple TTree with the readMuDst maker is described [here](/star-analysis/tree/main?tab=readme-ov-file#making-a-simple-root-file-for-analysis).

|Script name    | Description       |
|:--------------|:------------------|
|fwd_sim_ana1.C |Gen. vs. rec tracks|
|fwd_sim_ana2.C |MIP analysis       |


