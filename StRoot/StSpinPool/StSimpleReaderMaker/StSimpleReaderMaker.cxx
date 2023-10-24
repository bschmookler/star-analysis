//Class StSimpleReaderMaker

#include "StSimpleReaderMaker.h"

#include <iostream>

#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuFcsHit.h"
#include "StMuDSTMaker/COMMON/StMuFcsCollection.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"

#include "StFcsDbMaker/StFcsDbMaker.h"
#include "StFcsDbMaker/StFcsDb.h"

#include "StRoot/StEpdUtil/StEpdGeom.h"

#include "TTree.h"
#include "TFile.h"
#include "TObjArray.h"

ClassImp(StSimpleReaderMaker)                   // Macro for CINT compatibility

StSimpleReaderMaker::StSimpleReaderMaker( StMuDstMaker* maker ) : StMaker("StSimpleReaderMaker")
{ // Initialize and/or zero all public/private data members here.
  out_file = NULL;
  out_tree = NULL;
  mOutputFileName = "";
  mEventsProcessed = 0;

  mEpdgeo = new StEpdGeom;  //EPD geom.

  mMuDstMaker = maker ;     // Pass MuDst pointer to DstAnlysisMaker Class member functions
}

StSimpleReaderMaker::~StSimpleReaderMaker() 
{ // Destroy and/or zero out all public/private data members here.
}

Int_t StSimpleReaderMaker::Init( )
{ // Do once at the start of the analysis

  //FCS DB
  mFcsDb = static_cast<StFcsDb*>(GetDataSet("fcsDb"));
  //mFcsDb->setDbAccess(0);
  if (!mFcsDb) {
  	LOG_ERROR << "StSimpleReaderMaker::InitRun Failed to get StFcsDb" << endm;
        return kStFatal;
  }

  out_file = new TFile(mOutputFileName,"RECREATE");
  out_tree = new TTree("data","Simple Data Tree");

  out_tree->Branch("Cal_nhits",&Cal_nhits,"Cal_nhits/I");

  out_tree->Branch("Cal_detid",Cal_detid,"Cal_detid[Cal_nhits]/I");
  out_tree->Branch("Cal_hitid",Cal_hitid,"Cal_hitid[Cal_nhits]/I");
  out_tree->Branch("Cal_adcsum",Cal_adcsum,"Cal_adcsum[Cal_nhits]/I");
  out_tree->Branch("Cal_hit_energy",Cal_hit_energy,"Cal_hit_energy[Cal_nhits]/F");
  out_tree->Branch("Cal_hit_posx",Cal_hit_posx,"Cal_hit_posx[Cal_nhits]/F");
  out_tree->Branch("Cal_hit_posy",Cal_hit_posy,"Cal_hit_posy[Cal_nhits]/F");
  out_tree->Branch("Cal_hit_posz",Cal_hit_posz,"Cal_hit_posz[Cal_nhits]/F");

  return kStOK ; 

}

Int_t StSimpleReaderMaker::Make( ) 
{ // Do each event

  //Reset variables
  Cal_nhits = 0;

  // Get 'event' data 
  //StMuEvent* muEvent = mMuDstMaker->muDst()->event() ;

  // Get FCS data and fill TTree
  StMuFcsCollection* hits = mMuDstMaker->muDst()->muFcsCollection();  // Array containing the FCS Hits

  for(UInt_t ihit = 0 ; ihit < hits->numberOfHits() ; ihit++){

  	StMuFcsHit* hit = hits->getHit(ihit); // Pointer to a hit
	
	int det = hit->detectorId();

        if( det < kFcsNDet ){ //Some entries in MuDST file have detid > 5
    
    		Cal_detid[Cal_nhits] = det;
   		Cal_hitid[Cal_nhits] = hit->id();
    		Cal_adcsum[Cal_nhits] = hit->adcSum();
    		Cal_hit_energy[Cal_nhits] = hit->energy();

		if( det <= kFcsHcalSouthDetId ){
        		StThreeVectorD xyz = mFcsDb->getStarXYZ(det,hit->id());

    			Cal_hit_posx[Cal_nhits] = xyz.x();
    			Cal_hit_posy[Cal_nhits] = xyz.y();
    			Cal_hit_posz[Cal_nhits] = xyz.z();
		}
	      	else if(det==kFcsPresNorthDetId || det==kFcsPresSouthDetId){ //EPD as Pres.
			//Adapted from code StFcsEventDisplay.cxx
			double zepd=375.0;
			double zfcs=710.0+13.90+15.0;
			double zr=zfcs/zepd;
			int pp,tt,n;
			double x[5],y[5];
			double xsum(0), ysum(0);
	       	 	mFcsDb->getEPDfromId(det,hit->id(),pp,tt);
			mEpdgeo->GetCorners(100*pp+tt,&n,x,y);
 		
			//Get average of corner positions
			//N.B. Number of corners is usually 4, sometimes 5	
			for(int i=0; i<n; i++){
			    xsum += zr*x[i];
			    ysum += zr*y[i];
			}
			Cal_hit_posx[Cal_nhits] = xsum/n;
                	Cal_hit_posy[Cal_nhits] = ysum/n;
                	Cal_hit_posz[Cal_nhits] = zepd;
	    	}
		//Increment number of hits
		Cal_nhits++;
	}
  } //Loop over FCS hits

  mEventsProcessed++ ;
  
  //Fill TTree
  out_tree->Fill();
  return kStOK ;
  
}

Int_t StSimpleReaderMaker::Finish( )
{ // Do once at the end the analysis

  out_file->Write();
  out_file->Close();

  cout << "Total Events Processed in DstMaker " << mEventsProcessed << endl ;

  return kStOk ;  

}
