# Simulation with forward tracking

How to setup
------------
The forward tracking should work with the standard <i>stardev</i> environment. First build the geometry using the [build_geom.C](build_geom.C) file in this repository.

Single-particle simulation
--------------------------
To run a single-particle simulation, we can use the [runSimFlat.C](runSimFlat.C) and the [runSimBfc.C](runSimBfc.C) codes. To run 100 single negative pion events for pions with 30 GeV energy and Vz = 0, do the following:
```
root4star -b -q runSimFlat.C'(100,1,"pi-",30,0,0,1)'
```
This creates a ROOT file with the generated particles, as well as a <i>.fzd</i> file with the detector response information. 

The events can then be reconstructed by doing
```
root4star -b -q runSimBfc.C'(100,1,"pi-",202207,0,30)'
```
This generates a both MuDST file and a <i>.event</i> file. The MuDST file can then be processed to create a simple ROOT TTree for further analysis:
```
root4star -b -q 'readMudst.C(0,1,"input/pi-.MuDst.root")'
```

