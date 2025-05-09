# Analysis of production data
The production data with forward tracking is primarily being produced in PicoDst format. The scripts in this repository use PicoDst classes to perform data analysis.

How to setup
------------
Using the <i>stardev</i> enviroment, follow the instructions given [here](https://github.com/jdbrice/star-sw-1/wiki#accessing-up-to-date-code) to install the required PicoDst libraries. (Once the official PicoDst libraries are updated with the Forward classes, this step will not be needed.)

Creating a list of PicoDst files to analyze
--------------------------------------------
Run the following command to create a list of 60 PicoDst files from the test production run on April 30th, 2025:

```
ls -1 /gpfs01/star/pwg_tasks/FwdCalib/PROD/forwardCrossSection_2022/30042025/*.root | head -60 > infiles_30042025.lis
```

Analysis 1: MIP Peak analysis
------------------------------
This analysis code will use tracks projected to the front face of the FCS ECal and HCal to search for a minimum-ionizing particle (MIP) peak. Run as:

```
mkdir plots
root4star -l -b -q pico.C(2,"infiles_30042025.lis")
```


