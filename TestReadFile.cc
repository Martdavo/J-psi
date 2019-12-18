/* 
 * File:   TestReadFile.cc
 * Author: rafopar
 *
 * Created on September 23, 2019, 8:11 AM
 */


// ===== Hipo headers =====
#include <reader.h>
#include <dictionary.h>

// ===== Root headers =====
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>

// ====== C++ standard headers ======
#include <cstdlib>
#include <iostream>


using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

  const double Eb = 10.5;
  const double Mp = 0.9383;
  const double Mn = 0.93956;
  const double Md = 1.8756129;
  const double Me = 0.00051099892;
  const int DET_HTCC = 15; // HTCC is In the Rec::Cherenkov Bank. One should use Detector=15 for HTCC

  // ==== Create an output root file to store Root objects, e.g. histograms, trees, graphs etc
  TFile *file_out = new TFile("TestAnalyze.root", "Recreate");

  TH2D *h_th_P1 = new TH2D("h_th_P1", "", 200, 0., Eb, 200, 0., 65);
  TH2D *h_SampFrac_P1 = new TH2D("h_SampFrac_P1", "", 200, 0., Eb, 200, 0., 0.5);
  TH2D *h_SampFrac_P2 = new TH2D("h_SampFrac_P2", "", 200, 0., Eb, 200, 0., 0.5);
  TH1D *h_nphe_em1 = new TH1D("h_nphe_em1", "", 200, 0., 35.);

  TH2D *h_SampFrac_P11 = new TH2D("h_SampFrac_P11", "", 200, 0., Eb, 200, 0., 0.5);
  TH2D *h_SampFrac_P12 = new TH2D("h_SampFrac_P12", "", 200, 0., Eb, 200, 0., 0.5);
  TH2D *h_SampFrac_P13 = new TH2D("h_SampFrac_P13", "", 200, 0., Eb, 200, 0., 0.5);
  TH2D *h_SampFrac_P14 = new TH2D("h_SampFrac_P14", "", 200, 0., Eb, 200, 0., 0.5);
  TH2D *h_SampFrac_P15 = new TH2D("h_SampFrac_P15", "", 200, 0., Eb, 200, 0., 0.5);
  TH2D *h_SampFrac_P16 = new TH2D("h_SampFrac_P16", "", 200, 0., Eb, 200, 0., 0.5);

  TH2D *h_SampFrac_lu11 = new TH2D("h_SampFrac_lu11", "", 200, 0., 300, 200, 0., 0.5);
  TH2D *h_SampFrac_lv11 = new TH2D("h_SampFrac_lv11", "", 200, 0., 300, 200, 0., 0.5);
  TH2D *h_SampFrac_lw11 = new TH2D("h_SampFrac_lw11", "", 200, 0., 300, 200, 0., 0.5);

  TH2D *h_SampFrac_lu21 = new TH2D("h_SampFrac_lu21", "", 200, 0., 450, 200, 0., 0.5);
  TH2D *h_SampFrac_lv21 = new TH2D("h_SampFrac_lv21", "", 200, 0., 450, 200, 0., 0.5);
  TH2D *h_SampFrac_lw21 = new TH2D("h_SampFrac_lw21", "", 200, 0., 450, 200, 0., 0.5);

  TH1D *h_inv_mass  = new TH1D("h_inv_mass",  "", 300, 0, 4);
  TH1D *h_open_angle = new TH1D("h_open_angle", "", 400, -40, 400);
  TH2D *h_invmass_openangle = new TH2D("h_invmass_openangl", "", 200, 0., 4, 400, 0., 80);

        
  TH1D *h_inv_mass2  = new TH1D("h_inv_mass2",  "", 300, 0, 4);
  TH1D *h_open_angle2 = new TH1D("h_open_angle2", "", 400, -40, 400);
  TH2D *h_invmass_openangle2 = new TH2D("h_invmass_openangl2", "", 200, 0., 4, 400, 0., 80);

    
  TH1D *h_phi     = new TH1D("h_phi",     "", 400, 0, 400);
  TH1D *h_phi_el  = new TH1D("h_phi_el",  "", 400, 0, 400);
  TH1D *h_phi_ps  = new TH1D("h_phi_ps",  "", 400, 0, 400);

  TH1D *h_EPcale      = new TH1D("h_EPcale",  "", 200, -1., -1.);
  TH1D *h_EECine      = new TH1D("h_EEcine",  "", 200, -1., -1.);
  TH2D *h_ECin_Pcal_e = new TH2D("h_E_P_e", "", 200, 0., 0.2, 200, 0., 0.2);

  TH1D *h_EPcalps      = new TH1D("h_EPcalps",  "", 200, -1., -1.);
  TH1D *h_EECinps     = new TH1D("h_EEcinps",  "", 200, -1., -1.);
  TH2D *h_ECin_Pcal_ps = new TH2D("h_E_P_ps", "", 200, 0., 0.2, 200, 0., 0.2);
  TH2D *h_ECin_Pcal_ps2 = new TH2D("h_E_P_ps2", "", 200, 0., 0.2, 200, 0., 0.2);

  TH1D *h_ps_chi      = new TH1D("h_ps_chi",  "", 200, -6., 6.);
  TH1D *h_el_chi      = new TH1D("h_el_chi",  "", 200, -6., 6.);


  TH1D *h_aac0     = new TH1D("h_aacounter0",  "", 200, -1., 4.);
  TH1D *h_bbc0     = new TH1D("h_bbcounter0",  "", 200, -1., 4.);
  TH1D *h_aac     = new TH1D("h_aacounter",  "", 200, -1., 4.);
  TH1D *h_bbc     = new TH1D("h_bbcounter",  "", 200, -1., 4.);
    

  std::cout << " reading file example program (HIPO) " << __cplusplus << std::endl;

  if (argc < 2) {
    std::cout << " *** please provide a file name..." << std::endl;
    exit(0);
  }

  hipo::reader reader;
  char inputFile[256];
  for(int iFile=1; iFile<argc; iFile++){
    sprintf(inputFile, "%s", argv[iFile]);
    printf("Processing input data file ---> %s\n", inputFile);


    reader.open(inputFile);

    hipo::dictionary factory;

    reader.readDictionary(factory);

    factory.show();

    hipo::event event;
    int counter = 0;

    hipo::bank PART(factory.getSchema("REC::Particle"));
    hipo::bank bRecCalo(factory.getSchema("REC::Calorimeter"));
    hipo::bank bRecCC(factory.getSchema("REC::Cherenkov"));

        

    
    // Loop over events
    while (reader.next() == true) {
      reader.read(event);

      map<int, int> indPCal;
      map<int, int> indECin;
      map<int, int> indECout;
      map<int, int> indHTCC;

      event.getStructure(bRecCalo);
      int n_RecCalo = bRecCalo.getRows();

      for (int icalo = 0; icalo < n_RecCalo; icalo++) {

	int pindex = bRecCalo.getInt("pindex", icalo);
	int layer = bRecCalo.getInt("layer", icalo);

	// PCal, ECin, and ECout are part of the REC::Calorimeter Bank
	// they are distinguished by the variable "layer" 1 (PCal) 4 (ECin) and 7 (ECout)

	// here we will create the Part->ECindex relation
	if (layer == 1) {
	  indPCal[pindex] = icalo + 1;
	} else if (layer == 4) {
	  indECin[pindex] = icalo + 1;
	} else if (layer == 7) {
	  indECout[pindex] = icalo + 1;
	}

      }


      event.getStructure(bRecCC);
      int nCC = bRecCC.getRows();

      for (int iCC = 0; iCC < nCC; iCC++) {

	int pindex = bRecCC.getInt("pindex", iCC);
	int det = bRecCC.getInt("detector", iCC);

	// == Make sure the detector is HTCC
	if (det == DET_HTCC) {

	  // make the relation between Part->HTCCindex
	  indHTCC[pindex] = iCC + 1;
	}


      }


      event.getStructure(PART);
      //PART.show();
      int npart = PART.getRows();

      int n_em = 0;
      int n_ep = 0;
      int n_p = 0;
      int a = 0;
      int b = 0;
      int aa = 0.;
      int bb = 0.;
      int c = 0;
      int Sigma ; 
      double Electron_sector;
      double Positron_sector; 
      TLorentzVector L_electron, L_pozitron, L_inv, L_miss, L_empty, L_neutron , L_proton , L_beam, L_target ;
      /////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////
      for (int ipart = 0; ipart < npart; ipart++) {

	// Getting some variables from the Part bank

	int charge = PART.getInt("charge", ipart);  // Charge of the particle
	int stat = PART.getInt("status", ipart);    // status of the particle
	//int  EC_sector = bRecCalo.getInt("sector", ipart+1);
	// if (TMath::Abs(stat) >= 2000 && TMath::Abs(stat) <= 4000){
                    
	int pid = PART.getInt("pid", ipart);        // PID of the particle
	float px = PART.getFloat("px", ipart);      // X component of particl'es momentum 
	float py = PART.getFloat("py", ipart);      // Y component of particl'es momentum 
	float pz = PART.getFloat("pz", ipart);      // Z component of particl'es momentum
	float chi2 = PART.getFloat("chi2pid", ipart); // Chi2 of assigned PID

	float p = sqrt(px * px + py * py + pz * pz);
	double th = acos(pz / p) * TMath::RadToDeg();
	double phi = 30 + atan2(py/p,px/p)*TMath::RadToDeg() ;

	if (phi <   0) phi=phi+360;
	if (phi > 360) phi=phi-360;
	int sector = int(phi/60.) ;
	sector = sector + 1; 
	h_phi -> Fill(phi);
	double EPcal = 0;
	double lu_Pcal = 0 ;
	double lv_Pcal = 0 ;
	double lw_Pcal = 0 ; 
	double EECin = 0;
	double EECout = 0;

	const int Partconst = 1;
            
	// Getting energy depositions in PCal, ECin and ECout of the particle
	if (indPCal[ipart] > 0) {
	  EPcal = bRecCalo.getFloat("energy", indPCal[ipart] -Partconst );
	  lu_Pcal = bRecCalo.getFloat("lu", indPCal[ipart] - Partconst);
	  lv_Pcal = bRecCalo.getFloat("lv", indPCal[ipart] - Partconst);
	  lw_Pcal = bRecCalo.getFloat("lw", indPCal[ipart] - Partconst);
	}
	if (indECin[ipart] > 0) {
	  EECin = bRecCalo.getFloat("energy", indECin[ipart] - Partconst);
	}
	if (indECout[ipart] > 0) {
	  EECout = bRecCalo.getFloat("energy", indECout[ipart] - Partconst);
	}

	double Etot = EPcal + EECin + EECout;
	double EEC = EECin + EECout;

	// Getting the number of photoelectrons detected in HTCC from this particle
	float nphe = 0;
	if (indHTCC[ipart] > 0) {
	  nphe = bRecCC.getFloat("nphe", indHTCC[ipart] - Partconst);
	}

        
        
	// For electrons (PID == 11) detected in Forward Detector (2000 < status < 40000) fill some histograms
	if (pid == 11 && EPcal>0.06 && nphe > 2 && TMath::Abs(stat) >= 2000 && TMath::Abs(stat) <= 4000 ) {
	  h_th_P1->Fill(p, th);
	  //h_SampFrac_P1->Fill(p, Etot/p);
	  h_nphe_em1->Fill(nphe);
	  double E_el = sqrt(p*p+Me*Me);
	  L_electron.SetPxPyPzE(px,py,pz,E_el);
	  h_phi_el -> Fill(phi);
	  a=a+1;	  
	  Electron_sector = sector;
                
	  h_SampFrac_P11->Fill(p, Etot/p);
	  h_SampFrac_lu11->Fill(lu_Pcal, Etot/p);
	  h_SampFrac_lv11->Fill(lv_Pcal, Etot/p);
	  h_SampFrac_lw11->Fill(lw_Pcal, Etot/p);
                
                
	  if( lu_Pcal>15 && lv_Pcal>15 && lw_Pcal > 15 && TMath::Abs(chi2)>=-3 && TMath::Abs(chi2)<=3 && EECin>0 ){
	    h_EPcale -> Fill(EPcal);
	    h_EECine -> Fill(EECin);
	    if(EECin/p>=0.2-EPcal/p){
	      h_ECin_Pcal_e -> Fill(EPcal/p,EECin/p);
            
	      h_SampFrac_P1->Fill(p, Etot/p);
	      aa=aa+1;
	    }
	  }
        
        }

        
        
        
        if (pid == -11 && nphe > 2 && TMath::Abs(stat) >= 2000 && TMath::Abs(stat) <= 4000 ) {
	  //  h_th_P1->Fill(p, th);
	  //  h_SampFrac_P1->Fill(p, Etot/p);
	  //  h_nphe_em1->Fill(nphe);

	  h_ECin_Pcal_ps -> Fill(EPcal/p,EECin/p);
	  double EECincut = 0.2 - EPcal/p ; 
        
	  if( EECin>0)
            {
            
	      h_EPcalps -> Fill(EPcal);
	      h_EECinps -> Fill(EECin);
        
	      double E_ps = sqrt(p*p+Me*Me);
	      L_pozitron.SetPxPyPzE(px,py,pz,E_ps);
	      h_phi_ps -> Fill(phi);
	      Positron_sector = sector ;
	      h_SampFrac_lu21->Fill(lu_Pcal, Etot/p);
	      h_SampFrac_lv21->Fill(lv_Pcal, Etot/p);
	      h_SampFrac_lw21->Fill(lw_Pcal, Etot/p);
	      b=b+1;
	      h_ps_chi -> Fill(chi2);
	      Sigma = 3;
	      if(lu_Pcal>15 &&  lv_Pcal>15 && lw_Pcal > 15 && TMath::Abs(chi2)>=-Sigma && TMath::Abs(chi2)<=Sigma ){
		if(EECin/p>=0.2-EPcal/p){
		  h_ECin_Pcal_ps2 -> Fill(EPcal/p,EECin/p);
        
		  h_SampFrac_P2->Fill(p, Etot/p);
        
		  bb=bb+1;
		}
	      }
            }
        }

        

        if (pid == 2212  && TMath::Abs(stat) >= 2000 && TMath::Abs(stat) <= 4000) {
	  //  h_th_P1->Fill(p, th);
	  //  h_SampFrac_P1->Fill(p, Etot/p);
	  //  h_nphe_em1->Fill(nphe);
	  double E_p = sqrt(p*p+Mp*Mp);
	  L_proton.SetPxPyPzE(px,py,pz,E_p);
        
	  //  h_phi_ps -> Fill(phi);
	  // Positron_sector = sector ;
	  c=c+1;
        }

      } //}

      if (bb>0) h_bbc0 -> Fill(bb);
      if (aa>0) h_aac0 -> Fill(aa);

      if (aa>0 && bb>0){
        h_bbc -> Fill(bb);
        h_aac -> Fill(aa);}

      if( a>0 && b>0 ){
        //if( a==1 && b==1  && Electron_sector != Positron_sector){
        L_inv = L_electron + L_pozitron ; // making invaraint Lorentz vector 
        double  inv_mass= L_inv.Mag();   // making invariant mass
        double open_angle = L_electron.Angle(L_pozitron.Vect())* TMath::RadToDeg(); // angle between electron and positron vectors 
        h_inv_mass  -> Fill(inv_mass);
        h_open_angle->Fill(open_angle); 
        
      }

        
        
      if( aa>0 && bb>0 ){
        
        
        // if( aa==1 && bb==1  && Electron_sector != Positron_sector){
        L_inv = L_electron + L_pozitron ; // making invaraint Lorentz vector 
        double  inv_mass2= L_inv.Mag();   // making invariant mass
        double open_angle2 = L_electron.Angle(L_pozitron.Vect())* TMath::RadToDeg(); // angle between electron and positron vectors 
        h_inv_mass2  -> Fill(inv_mass2);
        h_open_angle2->Fill(open_angle2); 
        
      }

      counter++;
    } // file event counter
  } // iFile counter

        
    // Save all objects in the gDirectory to the output file    
  gDirectory->Write();
  return 0;
}
