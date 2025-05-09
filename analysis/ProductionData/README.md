# Analysis of production data
The Run22 p-p data including the Forward tracking is primarily being produced in PicoDst format. The scripts in this repository use PicoDst classes to perform data analysis on the produced data.

How to setup
------------
Using the <i>stardev</i> enviroment, follow the instructions given [here](https://github.com/jdbrice/star-sw-1/wiki#accessing-up-to-date-code) to install the required PicoDst libraries. (Once the official STAR PicoDst libraries are updated to incorporate the Forward detector classes, this step will no longer be needed.)

Creating a list of PicoDst files to analyze
--------------------------------------------
Run the following command to create a list of 60 PicoDst files from the test production run on April 30th, 2025:

```
ls -1 /gpfs01/star/pwg_tasks/FwdCalib/PROD/forwardCrossSection_2022/30042025/*.root | head -60 > infiles_30042025.lis
```

Analysis 1: MIP Peak analysis
------------------------------
This analysis code uses tracks projected to the FCS ECal and HCal to search for minimum-ionizing particle (MIP) peaks. Run as:

```
mkdir plots
root4star -l -b -q pico.C(2,"infiles_30042025.lis")
```


