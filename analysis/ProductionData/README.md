# Analysis of production data
The Run22 p-p data including the Forward tracking is being produced in PicoDst and MuDst formats. The scripts in this repository use PicoDst classes to perform data analysis on the produced PicoDst data files.

How to set up
-------------
The <i>stardev</i> environment now incorporates the Forward detector classes. So, no local compilation of any PicoDst libraries is required. 

If you wish to use a local version, follow the instructions given [here](https://github.com/jdbrice/star-sw-1/wiki#accessing-up-to-date-code) to install the required PicoDst libraries.

List of PicoDst files to analyze
---------------------------------
To use the runlist that was generated as part of the production that started in December 2025, create as symbolic link to that list:

```
ln -s /star/data14/GRID/NFS_FileList/production_pp500_2022_P25ib_st_physics_FwdTrack_picodst.txt production_pp500_2022.list
```

Analysis 1: MIP Peak analysis
------------------------------
This analysis code uses tracks projected to the FCS ECal and HCal to search for minimum-ionizing particle (MIP) peaks. Run as:

```
mkdir plots
root4star -l -b -q 'pico.C(2,"production_pp500_2022.list")'
```

