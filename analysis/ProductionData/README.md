# Analysis of production data
The Run22 p-p data including the Forward tracking is primarily being produced in PicoDst format. The scripts in this repository use PicoDst classes to perform data analysis on the produced data.

How to setup
------------
The <i>stardev</i> enviroment now incorporates the Forward detector classes. So, no local compilation of any PicoDst libraries is required. 

If you wish to use a local version, follow the instructions given [here](https://github.com/jdbrice/star-sw-1/wiki#accessing-up-to-date-code) to install the required PicoDst libraries.

Creating a list of PicoDst files to analyze
--------------------------------------------
Run the following command to create a list of 10 PicoDst files from the test production run on August 27th, 2025:

```
ls -1 /gpfs01/star/pwg_tasks/FwdCalib/PROD/forwardCrossSection_2022/27082025/*.root | head -10 | tee infiles_27082025.list
```

Analysis 1: MIP Peak analysis
------------------------------
This analysis code uses tracks projected to the FCS ECal and HCal to search for minimum-ionizing particle (MIP) peaks. Run as:

```
mkdir plots
root4star -l -b -q 'pico.C(2,"infiles_27082025.list")'
```


