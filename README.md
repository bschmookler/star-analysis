# star-analysis

Making a simple ROOT file for analysis
--------------------------------------
For a MuDst file containing STAR Forward Upgrade data, we have two codes which can read the MuDst file and create a simple ROOT TTree.

1. For all data containted in the MuDst file, read the MuDst classes directly. This approach can be found in the macro [readMudst.C](readMudst.C).
   
   To run the macro on a single file, the following command can be used:
   ```
   root4star -b -q 'readMudst.C(0,1,"input/zfa_prod/st_physics_23072003_raw_1000002.MuDst.root")'
   ```
   To run the macro on a [list of files](input/filelist.list), the following command can be used:
   ```
   root4star -b -q 'readMudst.C(0,5,"filelist.list")'
   ```
   
2. For the FCS data in the MuDst file, read the time-dependent ADC signal and then use the <i>StEvent</i> classes. This approach can be found in the macro [runMudst.C](runMudst.C).

   To run the macro on a single file, the following command can be used:
   ```
   root4star -b -q 'runMudst.C(0,1,"input/zfa_prod/st_physics_23072003_raw_1000002.MuDst.root")'
   ```

Simulation
----------
Information on running the full simulation is provided in the [simulation](simulation) subdirectory.

Analysis examples
-----------------
Some example analyses using the produced ROOT file can be found in the [analysis](analysis) subdirectory. In addition, analysis scripts using the PicoDst files can also be found in this subdirectory.
