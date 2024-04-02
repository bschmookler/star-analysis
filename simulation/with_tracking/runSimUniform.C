// macro to instantiate the Geant3 from within
// STAR  C++  framework and get the starsim prompt
// To use it do
//  root4star starsim.C

class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class StarGenEvent;
StarGenEvent   *event       = 0;

class StarPrimaryMaker;
StarPrimaryMaker *primary = 0;

class StarKinematics;
StarKinematics *kinematics = 0;

// Set magnetic field (reversed full)
Float_t field = -5.0;
//Float_t field = 0.0;

// ----------------------------------------------------------------------------
void geometry( TString tag, Bool_t agml=true ){
  TString cmd = "DETP GEOM "; cmd += tag; cmd += Form(" field=%f", field);
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> LoadGeometry(cmd);
  //  if ( agml ) command("gexec $STAR_LIB/libxgeometry.so");
}
// ----------------------------------------------------------------------------
void command( TString cmd ){
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> Do( cmd );
}
// ----------------------------------------------------------------------------
void trig( Int_t n=1, char* pid="ele", int npart=1, int print=0){
    char* PID = pid;
    char* ele="e-";
    char* pos="e+";
    char *piplus="pi+";
    if(pid[0]=='e' && pid[1]=='l' && pid[2]=='e') PID=ele;
    if(pid[0]=='p' && pid[1]=='o' && pid[2]=='s') PID=pos;
    if(pid[0]=='p' && pid[1]=='i' && pid[2]=='p') PID=piplus;
    for ( Int_t i=0; i<n; i++ ) {
	double pi=3.141592654;
	cout << "==== Event="<<i<<"===="<<endl;
	// Clear the chain from the previous event
	chain->Clear();

	//Generate from E = [1,100]GeV and eta = [2.25,4.75]
	kinematics->SetAttr("energy",1);
	kinematics->Kine(npart, PID, 1., 100., 2.25, 4.75, 0.0, 2*pi);	    
	
	// Generate the event
	chain->Make();
	// Print the event
	if(print>0) primary->event()->Print();
	if(print>1) {
	    TIter Iterator = primary->event()->IterAll();
	    StarGenParticle *p = 0;
	    while( ( p = (StarGenParticle*)Iterator.Next() ) ){
	      TLorentzVector v = p->momentum();
	      cout << Form(" ===> e=%6.3f pt=%7.3f eta=%7.3f phi=%7.3f\n",v.E(),v.Pt(),v.Eta(),v.Phi());
	    }
	}
	if(print>2) command("gpri hits");
    }
}
// ----------------------------------------------------------------------------
void Kinematics(){
  gSystem->Load( "libKinematics.so");
  kinematics = new StarKinematics();    

  primary->AddGenerator(kinematics);
}
// ----------------------------------------------------------------------------
void runSimUniform( Int_t nevents=1000, Int_t run=1, char* pid="ele", int npart=1, int ecal=1, int print=0){ 

  gSystem->Load( "libStarRoot.so" );
  gROOT->SetMacroPath(".:/star-sw/StRoot/macros/:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");

  gROOT->ProcessLine(".L bfc.C");
  {
    TString simple = "sdt20211016 y2023 geant gstar agml usexgeom";
    bfc(0, simple );
  }

  gSystem->Load( "libVMC.so");
  gSystem->Load( "StarGeneratorUtil.so" );
  gSystem->Load( "StarGeneratorEvent.so" );
  gSystem->Load( "StarGeneratorBase.so" );
  gSystem->Load( "libMathMore.so"   );  
  gSystem->Load( "xgeometry.so"     );
  
  // Create the primary event generator and insert it
  // before the geant maker
  primary = new StarPrimaryMaker();
  {
    chain   -> AddBefore( "geant", primary );
    // Set the output filename for the event record
    primary -> SetFileName(Form("%s.run%i.root",pid,run));

    // Set the x,y,z vertex and distribution
    primary -> SetVertex( 0., 0., 0.);
    primary -> SetSigma ( 0., 0., 0. );
  }

  // Initialize primary event generator and all sub makers
  Kinematics();
  primary -> Init();

  // Initialize random number generator
  StarRandom::capture();
  StarRandom::seed(run);

  // Setup geometry and set starsim to use agusread for input
  //geometry("dev2022=1");geometry("dev2022");geometry("y2023");

  command("gkine -4 0");
  command(Form("gfile o %s.run%d.fzd",pid,run));

  if(ecal==0){
      cout << "TURNING OFF ECAL!!!"<<endl;
      command("DETP WCAL wver.active=0");
  }

  // Switch off some physics
  const Char_t *cmds[] = { 
    "CKOV 0"
  };
  for ( UInt_t i=0;i<sizeof(cmds)/sizeof(Char_t*); i++ ){command( cmds[i] );}
    
  // Trigger on nevents
  trig(nevents, pid, npart, print);

  command("call agexit");
}
// ----------------------------------------------------------------------------

