# SimpleTree: Making a Simple ROOT File for Analysis

## Setup

Compile the maker classes in [StRoot/StSpinPool](StRoot/StSpinPool):

```
cd StRoot
cons
```

---

## Running the Macros

For a MuDst file containing STAR Forward Upgrade data, two macros are available for reading the MuDst file and creating a simple ROOT TTree.

### 1. 📄 `readMudst.C` — Full MuDst data

Reads the MuDst classes directly, capturing all data contained in the MuDst file.

**Single file:**
```
root4star -b -q 'readMudst.C(0,1,"input/zfa_prod/st_physics_23072003_raw_1000002.MuDst.root")'
```

**List of files** (see [`input/filelist.list`](input/filelist.list)):
```
root4star -b -q 'readMudst.C(0,5,"filelist.list")'
```

### 2. 📄 `runMudst.C` — FCS data via StEvent

Reads the time-dependent ADC signal from the FCS data in the MuDst file, then uses the *StEvent* classes.

**Single file:**
```
root4star -b -q 'runMudst.C(0,1,"input/zfa_prod/st_physics_23072003_raw_1000002.MuDst.root")'
```

---

## Analysis Examples

Example analyses using the produced simple ROOT file can be found in the [analysis](../analysis) directory.
