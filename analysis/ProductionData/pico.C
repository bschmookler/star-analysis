void pico(int do=2, const Char_t *inFile = "infiles.lis") {
  gSystem->Load("libPhysics");
  gSystem->Load("St_base");
  gSystem->Load("StChain");
  gSystem->Load("St_Tables");
  gSystem->Load("StUtilities");      
  //gSystem->Load("StTreeMaker");
  //gSystem->Load("StIOMaker");
  //gSystem->Load("StBichsel");
  gSystem->Load("StEvent");
  //gSystem->Load("StTpcDb");
  //gSystem->Load("StEventUtilities");
  //gSystem->Load("StTriggerDataMaker");
  gSystem->Load("StDbLib");
  gSystem->Load("StFcsDbMaker");
  gSystem->Load("StPicoEvent");
  
  if(do==0) gROOT->ProcessLine(Form(".x match.C+(\"%s\")",inFile));
  if(do==1) gROOT->ProcessLine(Form(".x dilepton.C+(\"%s\")",inFile));
  if(do==2) gROOT->ProcessLine(Form(".x mip_ana.C+(\"%s\")",inFile));
}

