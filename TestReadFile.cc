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

  std::cout << " reading file example program (HIPO) " << __cplusplus << std::endl;

  char inputFile[256];

  if (argc > 1) {
    sprintf(inputFile, "%s", argv[1]);
    //sprintf(outputFile,"%s",argv[2]);
  } else {
    std::cout << " *** please provide a file name..." << std::endl;
    exit(0);
  }

  const double Eb = 10.5;
  const double Mp = 0.9383;
  const double Mn = 0.93956;
  const double Md = 1.8756129;
  const double Me = 0.00051099892;
  const int DET_HTCC = 15; // HTCC is In the Rec::Cherenkov Bank. One should use Detector=15 for HTCC

  hipo::reader reader;
  reader.open(inputFile);

  hipo::dictionary factory;

  reader.readDictionary(factory);

  factory.show();

  hipo::event event;
  int counter = 0;

  hipo::bank PART(factory.getSchema("REC::Particle"));
  hipo::bank bRecCalo(factory.getSchema("REC::Calorimeter"));
  hipo::bank bRecCC(factory.getSchema("REC::Cherenkov"));

    
  // ==== Create an output root file to store Root objects, e.g. histograms, trees, graphs etc
  TFile *file_out = new TFile("TestAnalyze.root", "Recreate");

  TH2D *h_th_P1 = new TH2D("h_th_P1", "", 200, 0., Eb, 200, 0., 65);
  TH2D *h_SampFrac_P1 = new TH2D("h_SampFrac_P1", "", 200, 0., Eb, 200, 0., 0.5);
  TH1D *h_nphe_em1 = new TH1D("h_nphe_em1", "", 200, 0., 35.);

  TH1D *h_inv_mass  = new TH1D("h_inv_mass",  "", 300, 0, 4);
  TH1D *h_open_angle = new TH1D("h_open_angle", "", 400, -40, 400);
  TH2D *h_invmass_openangle = new TH2D("h_invmass_openangl", "", 200, 0., 4, 400, 0., 80);
  

 
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
    double Electron_sector;
    double Positron_sector; 
    TLorentzVector L_electron, L_pozitron, L_inv, L_miss, L_empty, L_neutron , L_proton , L_beam, L_target ;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (int ipart = 0; ipart < npart; ipart++) {

      // Getting some variables from the Part bank

      int charge = PART.getInt("charge", ipart);  // Charge of the particle
      int stat = PART.getInt("status", ipart);    // status of the particle
      if (TMath::Abs(stat) >= 2000 && TMath::Abs(stat) <= 4000){
            
	int pid = PART.getInt("pid", ipart);        // PID of the particle
	float px = PART.getFloat("px", ipart);      // X component of particl'es momentum 
	float py = PART.getFloat("py", ipart);      // Y component of particl'es momentum 
	float pz = PART.getFloat("pz", ipart);      // Z component of particl'es momentum      

      	float p = sqrt(px * px + py * py + pz * pz);
	double th = acos(pz / p) * TMath::RadToDeg();
	double phi = atan2(py/p,px/p) * TMath::RadToDeg(); 
	if(phi <   0){phi=phi+360;}
	if(phi > 360){phi=phi-360;}
	int sector = phi/60. ;
	sector = sector + 1; 

	double EPcal = 0;
	double EECin = 0;
	double EECout = 0;
	
	// Getting energy depositions in PCal, ECin and ECout of the particle
	if (indPCal[ipart] > 0) {
	  EPcal = bRecCalo.getFloat("energy", indPCal[ipart] - 1);
	}
	if (indECin[ipart] > 0) {
	  EECin = bRecCalo.getFloat("energy", indECin[ipart] - 1);
	}
	if (indECout[ipart] > 0) {
	  EECout = bRecCalo.getFloat("energy", indECout[ipart] - 1);
	}

	double Etot = EPcal + EECin + EECout;
	double EEC = EECin + EECout;

	// Getting the number of photoelectrons detected in HTCC from this particle
	float nphe = 0;
	if (indHTCC[ipart] > 0) {
	  nphe = bRecCC.getFloat("nphe", indHTCC[ipart] - 1);
	}
            
	// For electrons (PID == 11) detected in Forward Detector (2000 < status < 40000) fill some histograms
	if (pid == 11 && TMath::Abs(stat) >= 2000 && TMath::Abs(stat) <= 4000) {
	  h_th_P1->Fill(p, th);
	  h_SampFrac_P1->Fill(p, Etot/p);
	  h_nphe_em1->Fill(nphe);
	  double E_el = sqrt(p*p+Me*Me);
	  L_electron.SetPxPyPzE(px,py,pz,E_el);
	  Electron_sector = sector; 
	  a=a+1;
	  
	}

	if (pid == -11 && TMath::Abs(stat) >= 2000 && TMath::Abs(stat) <= 4000) {
	  h_th_P1->Fill(p, th);
	  h_SampFrac_P1->Fill(p, Etot/p);
	  h_nphe_em1->Fill(nphe);
	  double E_ps = sqrt(p*p+Me*Me);
	  L_pozitron.SetPxPyPzE(px,py,pz,E_ps);
	  Positron_sector = sector ;
	  b=b+1;
	}

      } }

    if( a==1 && b==1 && Electron_sector != Positron_sector){
      L_inv = L_electron + L_pozitron ; // making invaraint Lorentz vector 
      double  inv_mass= L_inv.Mag2();   // making invariant mass
      double open_angle = L_electron.Angle(L_pozitron.Vect())* TMath::RadToDeg(); // angle between electron and positron vectors 
      h_inv_mass  -> Fill(inv_mass);
      h_open_angle->Fill(open_angle); 
      
    }

    counter++;
  }

    
  // Save all objects in the gDirectory to the output file    
  gDirectory->Write();
  return 0;
}
