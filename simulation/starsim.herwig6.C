// macro to instantiate the Geant3 from within
// STAR  C++  framework and get the starsim prompt
// To use it do
//  root4star starsim.C

class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class StarGenEvent;
StarGenEvent   *event       = 0;

class StarPrimaryMaker;
StarPrimaryMaker *_primary   = 0;

//class StarHerwig6;
//StarHerwig6 *gHerwig        = 0;

// ----------------------------------------------------------------------------
void geometry( TString tag, Bool_t agml=true )
{
  TString cmd = "DETP GEOM "; cmd += tag;
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> LoadGeometry(cmd);
  //  if ( agml ) command("gexec $STAR_LIB/libxgeometry.so");
}
// ----------------------------------------------------------------------------
void command( TString cmd )
{
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> Do( cmd );
}
// ----------------------------------------------------------------------------
void trig( Int_t n=1 )
{
  for ( Int_t i=0; i<n; i++ ) {
    chain->Clear();
    chain->Make();
    _primary->event()->Print();     // Print the primary event
  }
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void Herwig6( TString mode="pp" )
{
  
  gSystem->Load( "libHerwig6_5_20.so");

  //StarHerwig6 *herwig6 = new StarHerwig6("herwig6");
  herwig6 = new StarHerwig6("herwig6");
  if ( mode == "pp" )
  {
    Double_t pblue[]={0.,0.,250.0};
    Double_t pyell[]={0.,0.,-250.0};
    herwig6->SetFrame("3MOM", pblue, pyell );
    herwig6->SetBlue("proton");
    herwig6->SetYell("proton");
    herwig6->SetProcess(1000);
  }
    
  _primary->AddGenerator(herwig6);
}
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
void starsim( Int_t nevents=1, Int_t rngSeed=1234 )
{ 

  gROOT->ProcessLine(".L bfc.C");
  {
    TString simple = "y2023 geant gstar usexgeom agml ";
    bfc(0, simple );
  }

  gSystem->Load( "libVMC.so");

  gSystem->Load( "StarGeneratorUtil.so" );
  gSystem->Load( "StarGeneratorEvent.so" );
  gSystem->Load( "StarGeneratorBase.so" );

  gSystem->Load( "libMathMore.so"   );  
  gSystem->Load( "xgeometry.so"     );

  // Setup RNG seed and captuire ROOT TRandom
  StarRandom::seed(rngSeed);
  StarRandom::capture();
  
  //
  // Create the primary event generator and insert it
  // before the geant maker
  //
  //  StarPrimaryMaker *
  _primary = new StarPrimaryMaker();
  {
    _primary -> SetFileName( "herwig6.starsim.root");
    _primary -> SetVertex( 0.1, -0.1, 0.0 );
    _primary -> SetSigma ( 0.1,  0.1, 30.0 );
    chain -> AddBefore( "geant", _primary );
  }

  //
  // Setup an event generator
  //
  Herwig6( "pp" );

  //
  // Setup cuts on which particles get passed to geant for
  // simulation.  (To run generator in standalone mode,
  // set ptmin=1.0E9.)
  //                     ptmin  ptmax
   _primary->SetPtRange  (0.0,  -1.0);         // GeV
  //                     etamin etamax
  _primary->SetEtaRange ( +1.0, +5.0 );
  //                     phimin phimax
  _primary->SetPhiRange ( 0., TMath::TwoPi() );

  // 
  // Setup a realistic z-vertex distribution:
  //   x = 0 gauss width = 1mm
  //   y = 0 gauss width = 1mm
  //   z = 0 gauss width = 30cm
  //
  _primary->SetVertex( 0., 0., 0. );
  _primary->SetSigma( 0.1, 0.1, 30.0 );


  //
  // Initialize primary event generator and all sub makers
  //
  _primary -> Init();

  //
  // Setup geometry and set starsim to use agusread for input
  //
  geometry("y2023");
  command("gkine -4 0");
  command("gfile o herwig6.starsim.fzd");
  

  //
  // Trigger on nevents
  //
  trig( nevents );

  command("call agexit");  // Make sure that STARSIM exits properly

}
// ----------------------------------------------------------------------------
