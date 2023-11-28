# STAR Simulation

How to setup
------------
Using the <i>stardev</i> environment, follow the instructions on [this site](https://www.star.bnl.gov/protected/spin/akio/fcs/howto_MC_github.html). (N.B. I decided to use the main branch to download the files, instead of the other one mentioned in the link.)

In order to save the simulation output to a MuDST file, some changes need to be made to the runSimBfc.C file. These changes are implemented in the [runSimBfc.C](runSimBfc.C) file in this repository.

For the FCS, we extract the hit energy directly from the <i>.fzd</i> file by running the <i>WaveFormFitMaker</i> with the option given [here](https://github.com/star-bnl/star-sw/blob/main/StRoot/StFcsWaveformFitMaker/StFcsWaveformFitMaker.cxx#L475). In order to save the ADC sum to the output MuDST file, we need to add the line ```hits[i]->setAdcSum(adc);``` right after the code [here](https://github.com/star-bnl/star-sw/blob/main/StRoot/StFcsFastSimulatorMaker/StFcsFastSimulatorMaker.cxx#L323).

It is also possible to run using the libraries already in the <i>stardev</i> environment, provided you still compile the libraries in ```StSpinPool```. If you wish to save the ADC sum you need to modify the ```StFcsFastSimulatorMaker``` as described above.

Single-particle simulation
--------------------------
To run a single-particle simulation, we can use the [runSimFlat.C](runSimFlat.C) and the [runSimBfc.C](runSimBfc.C) codes. To run 100 single negative pion events for pions with 30 GeV energy and Vz = 0, do the following:
```
root4star -b -q runSimFlat.C'(100,1,"pi-",30,0,0,1)'
```
Note how the code contains the following lines:
```
if(e>0.0)
{
	kinematics->SetAttr("energy",1);
	kinematics->Kine(npart, PID, e-0.01, e+0.01, 3.0,  3.01, 0.0, 2*pi);
}
```
This means the pion will be generated within the pseudorapidity range of 3.0 < eta < 3.01.

This creates a ROOT file with the generated particles, as well as a <i>.fzd</i> file with the detector response information. The events can then be reconstructed by doing
```
root4star -b -q runSimBfc.C'(100,1,"pi-",202207,0,30)'
```
This generates a MuDST file, which can then be processed to create a simple ROOT TTree for further analysis:
```
root4star -b -q 'readMudst.C(0,1,"input/pi-.MuDst.root")'
```

### Single-particle simulation with tracking
Follow the instructions in [this](with_tracking) subdirectory in order to include the forward tracking in the reconstruction.

Pythia8 simulation
------------------
To generate events and run those events through the STAR detector simulation, we use the [starsim.pythia8.C](starsim.pythia8.C) code. We create a softlink called ```starsim.C``` and run as follows:
```
root4star -b -q 'starsim.C(1000)'
```
This creates a ROOT file with the generated particles, as well as a <i>.fzd</i> file with the detector response information. The events can then be reconstructed by doing
```
root4star -b -q 'runSimBfc.C(1000,1,"pythia8")'
```
This generates a MuDST file, which can then be processed to create a simple ROOT TTree for further analysis:
```
root4star -b -q 'readMudst.C(0,1,"input/pythia8.MuDst.root")'
```

### Pythia8 simulation with event filtering
You can also filter the generated Pythia8 events before passing them through the detector simulation. In order to include the [FCS jet filter](https://github.com/star-bnl/star-sw/blob/main/StRoot/StarGenerator/FILT/FcsJetFilter.cxx), use the [starsim.pythia8_filter.C](starsim.pythia8_filter.C) code in this repository. All other steps are the same as described above.

Pythia6 simulation
------------------
To use the Pythia6 event generator instead of Pythia8, have the ```starsim.C``` softlink point to the [starsim.pythia6.C](starsim.pythia6.C) code. Then run as follows:
```
root4star -b -q 'starsim.C(1000)'
root4star -b -q 'runSimBfc.C(1000,1,"pythia6")'
root4star -b -q 'readMudst.C(0,1,"input/pythia6.MuDst.root")'
```

Herwig6 simulation
------------------
To use the Herwig6 event generator instead of Pythia8, have the ```starsim.C``` softlink point to the [starsim.herwig6.C](starsim.herwig6.C) code. Then run as follows:
```
root4star -b -q 'starsim.C(1000)'
root4star -b -q 'runSimBfc.C(1000,1,"herwig")'
root4star -b -q 'readMudst.C(0,1,"input/herwig6.MuDst.root")'
```
<strong> For Herwig6, the particle kinematics appear to be generated incorrectly. Need to look into this at some point. </strong>

Running higher statistics simulations using the batch farm
-----------------------------------------------------------
To run events using the STAR scheduler, follow the instructions in [this](job_submission) subdirectory.
