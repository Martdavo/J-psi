// ===CVT sos program====

// ===== Hipo headers =====
#include <reader.h>
#include <dictionary.h>
//#include <header.h>
//#include <event.h>

// deutron pid == 45

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

const int Partconst = 1 ; 
const double Mp = 0.9383; //for RGA
const double Mn = 0.93956;
const double Md = 1.8756129; //for RGB
const double Me = 0.00051099892;
const double Mpion  = 0.13957;
const double Mpion0 = 0.13497;
const double Mtarget = Mp;
const int DET_HTCC = 15; // HTCC is In the Rec::Cherenkov Bank. One should use Detector=15 for HTCC
const int DET_FT = 10;
const double Eb =10.6;
const double charge_e  =1.60217662*pow(10,-10);
TLorentzVector L_beam, L_target ;


double gen_inv_mass  (TLorentzVector L_inv_pt1,TLorentzVector L_inv_pt2);
double gen_trmom_x   (TLorentzVector L_inv_pt1,TLorentzVector L_inv_pt2,TLorentzVector L_inv_pt3);
double gen_trmom_y   (TLorentzVector L_inv_pt1,TLorentzVector L_inv_pt2,TLorentzVector L_inv_pt3);
double gen_trmom_z   (TLorentzVector L_inv_pt1,TLorentzVector L_inv_pt2,TLorentzVector L_inv_pt3);
double gen_trmom     (TLorentzVector L_inv_pt1,TLorentzVector L_inv_pt2,TLorentzVector L_inv_pt3);
double P_ee          (double px1,double py1,double pz1,double px2,double py2,double pz2);
double gen_miss_mass (double m, TLorentzVector L_miss_pt1,TLorentzVector L_miss_pt2,TLorentzVector L_miss_pt3);
double gen_miss_mass0(double m, TLorentzVector L_miss_pt1,TLorentzVector L_miss_pt2,TLorentzVector L_miss_pt3);
double gen_gamma_E   (double m, TLorentzVector L_miss_pt1,TLorentzVector L_miss_pt2,TLorentzVector L_miss_pt3);


 
int main(int argc, char** argv) {
  ofstream out_data("out_data_CVTsosRGK2_5693.dat");
  // ==== Create an output root file to store Root objects, e.g. histograms, trees, graphs etc
  TFile *file_out = new TFile("CVT_sosRGK2_5693.root", "Recreate");

  // -- histograms --
  TH1D *h_Electron_number_s1    = new TH1D("h_Electron_number_s1",  "", 200, 0, 10.);
  TH1D *h_Electron_number_s2    = new TH1D("h_Electron_number_s2",  "", 200, 0, 10.);
  TH1D *h_Electron_number_s3    = new TH1D("h_Electron_number_s3",  "", 200, 0, 10.);
  TH1D *h_Positron_number       = new TH1D("h_Positron_number",     "", 200, 0, 10.);
  TH1D *h_Positive_CD_number    = new TH1D("h_Positive_CD_number",  "", 200, 0, 10.);
  TH1D *h_Negative_CD_number    = new TH1D("h_Negative_CD_number",  "", 200, 0, 10.);

  TH1D *h_Electron_th       = new TH1D("h_Electron_th",     "", 200, 0, 100.);
  TH1D *h_Positron_th       = new TH1D("h_Positron_th",     "", 200, 0, 100.);
  TH1D *h_Positive_CD_th    = new TH1D("h_Positive_CD_th",  "", 200, 0, 100.);
  TH1D *h_Negative_CD_th    = new TH1D("h_Negative_CD_th",  "", 200, 0, 100.);

  TH1D *h_Positive_CD_th1    = new TH1D("h_Positive_CD_th1",  "", 200, 0, 100.);
  TH1D *h_Negative_CD_th1    = new TH1D("h_Negative_CD_th1",  "", 200, 0, 100.);
  TH1D *h_Positive_CD_th2    = new TH1D("h_Positive_CD_th2",  "", 200, 0, 100.);
  TH1D *h_Negative_CD_th2    = new TH1D("h_Negative_CD_th2",  "", 200, 0, 100.);
  TH2D *h_positive_angles    = new TH2D("h_positive_angles",  "", 200, 0, 180, 200, 0, 360.);
  TH2D *h_negative_angles    = new TH2D("h_negative_angles",  "", 200, 0, 180, 200, 0.,360.);


  
  TH2D *h_bullet_building            = new TH2D("h_bullet_building",        "", 200, 0, 100, 200, 0, 1.);
  TH2D *h_positive_angle_dp          = new TH2D("h_positive_angle_dp",         "", 200, -180, 180, 200, -180, 180.);
  TH2D *h_negative_angle_dp          = new TH2D("h_negative_angle_dp",         "", 200, -180, 180, 200, -180, 180.);
  TH2D *h_positive_angle_dp_s1       = new TH2D("h_positive_angle_dp_s1",      "", 200, -180, 180, 200, -180, 180.);
  TH2D *h_negative_angle_dp_s1       = new TH2D("h_negative_angle_dp_s1",      "", 200, -180, 180, 200, -180, 180.);
  TH2D *h_positive_angle_dp_s2       = new TH2D("h_positive_angle_dp_s2",      "", 200, -180, 180, 200, -180, 180.);
  TH2D *h_negative_angle_dp_s2       = new TH2D("h_negative_angle_dp_s2",      "", 200, -180, 180, 200, -180, 180.);
  TH2D *h_positive_angle_dp_s3       = new TH2D("h_positive_angle_dp_s3",      "", 200, -180, 180, 200, -180, 180.);
  TH2D *h_negative_angle_dp_s3       = new TH2D("h_negative_angle_dp_s3",      "", 200, -180, 180, 200, -180, 180.);

  TH2D *h_positive_angle_dp_cut   = new TH2D("h_positive_angle_dp_cut",  "", 200, -10, 10, 200, -10, 10.);
  TH2D *h_negative_angle_dp_cut   = new TH2D("h_negative_angle_dp_cut",  "", 200, -10, 10, 200, -10, 10.);
  
  TH2D *h_positive_angle_dp_cut_s1   = new TH2D("h_positive_angle_dp_cut_s1",  "", 200, -10, 10, 200, -10, 10.);  
  TH2D *h_negative_angle_dp_cut_s1   = new TH2D("h_negative_angle_dp_cut_s1",  "", 200, -10, 10, 200, -10, 10.);
  TH2D *h_positive_angle_dp_cut_s2   = new TH2D("h_positive_angle_dp_cut_s2",  "", 200, -10, 10, 200, -10, 10.);  
  TH2D *h_negative_angle_dp_cut_s2   = new TH2D("h_negative_angle_dp_cut_s2",  "", 200, -10, 10, 200, -10, 10.);
  TH2D *h_positive_angle_dp_cut_s3   = new TH2D("h_positive_angle_dp_cut_s3",  "", 200, -10, 10, 200, -10, 10.);  
  TH2D *h_negative_angle_dp_cut_s3   = new TH2D("h_negative_angle_dp_cut_s3",  "", 200, -10, 10, 200, -10, 10.);
  

  TH1D *h_FitMed      = new TH1D("h_FitMed",      "", 400, -2, 10.);
  TH1D *h_FitMed2     = new TH1D("h_FitMed2",     "", 400, -2, 10.);  
  TH1D *h_FitMed3     = new TH1D("h_FitMed3",     "", 400, -2, 10.);


  TH1D *h_BMT_ID          = new TH1D("h_BMT_ID0",      "", 400, 500., 2000.);
  TH1D *h_BMT_sector      = new TH1D("h_BMT_sector",   "", 4,   -0.5, 3.5);
  TH1D *h_BMT_Region      = new TH1D("h_BMT_Region",   "", 10,  0.,   10.);
  TH1D *h_BMT_Layer       = new TH1D("h_BMT_Layer",    "", 10,  0.,   10.);
  
  TH1D *h_chi2_pos    = new TH1D("h_chi2_pos",  "", 400, -20, 20.);
  TH1D *h_chi2_neg    = new TH1D("h_chi2_neg",  "", 400, -20, 20.);

  // -end of histograms --

  int el_counter=0,el_counter_s1=0,el_counter_s2=0,el_counter_s3=0;
  int pos_counter=0, pos_counter_s1=0,pos_counter_s2=0,pos_counter_s3=0;
  int neg_counter=0, neg_counter_s1=0,neg_counter_s2=0,neg_counter_s3=0;
  int pos_counter_cut=0, pos_counter_cut_s1=0,pos_counter_cut_s2=0,pos_counter_cut_s3=0;
  int neg_counter_cut=0, neg_counter_cut_s1=0,neg_counter_cut_s2=0,neg_counter_cut_s3=0;
  int N_run=0, N_event=0;
  int event_counter=0;

  int n0=0,n1=0,n2=0,n3=0,n4=0,n5=0,n6=0,n7=0,n8=0,n9=0,nab=0;

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
  
    // hipo::header header;

    hipo::bank PART    (factory.getSchema("REC::Particle"));
    hipo::bank bRecCalo(factory.getSchema("REC::Calorimeter"));
    hipo::bank bRecCC  (factory.getSchema("REC::Cherenkov"));
    hipo::bank bEvent  (factory.getSchema("REC::Event"));
    hipo::bank bconfig (factory.getSchema("RUN::config"));
    hipo::bank bFowTgg (factory.getSchema("REC::ForwardTagger"));
    hipo::bank bscaler (factory.getSchema("RUN::scaler"));
    hipo::bank bRecSC  (factory.getSchema("REC::Scintillator"));
    hipo::bank bBmtCr  (factory.getSchema("BMTRec::Crosses"));
    hipo::bank bcvtTr  (factory.getSchema("CVTRec::Tracks"));
        
    
    // Loop over events
    while (reader.next() == true) {
      reader.read(event);

      map<int, int> indPCal;
      map<int, int> indECin;
      map<int, int> indECout;
      map<int, int> indHTCC;
      map<int, int> indRCON;
      map<int, int> indFT;
      map<int, int> indSC;
      map<int, int> indCVT;
      map<int, int> indBMT; 
      
      event_counter++;
      event.getStructure(bEvent);
      int npart = PART.getRows();
      int Pid;

      int CVT_FM1=0,CVT_FM2=0;
      int Electron_number=0, Positron_number=0, Negative_CD_number=0,Positive_CD_number=0, Negative_FD_number=0,Positive_FD_number=0 ;
      int Electron_sector[300],Positron_sector[300],Negative_CD_sector[300],Positive_CD_sector[300],Negative_FD_sector[300],Positive_FD_sector[300];;
      //TLorentzVector L_electron[300],L_positive_CD[300];  
      double P_el[300],P_elx[300],P_ely[300],P_elz[300],vt_el[300],el_phi[300],el_theta[300];
      double P_ps[300],P_psx[300],P_psy[300],P_psz[300],vt_ps[300],ps_phi[300],ps_theta[300];
      double P_CD_pos[300],P_CD_posx[300],P_CD_posy[300],P_CD_posz[300],vt_CD_pos[300], pos_CD_phi[300],pos_CD_theta[300];
      double P_FD_pos[300],P_FD_posx[300],P_FD_posy[300],P_FD_posz[300],vt_FD_pos[300], pos_FD_phi[300],pos_FD_theta[300];
      double P_CD_neg[300],P_CD_negx[300],P_CD_negy[300],P_CD_negz[300],vt_CD_neg[300], neg_CD_phi[300],neg_CD_theta[300];
      double P_FD_neg[300],P_FD_negx[300],P_FD_negy[300],P_FD_negz[300],vt_FD_neg[300], neg_FD_phi[300],neg_FD_theta[300];
      
      double CFC = bEvent.getFloat("beamCharge",0);
     
      event.getStructure(bconfig);
      N_event = bconfig.getInt("event",0);
      N_run   = bconfig.getInt("run",0);

      event.getStructure(bscaler);
      double fcupg = bscaler.getFloat("fcupgated",0);
      double fcup  = bscaler.getFloat("fcup",0);
      
      event.getStructure(bRecCalo);
      int n_RecCalo = bRecCalo.getRows();

      
      

      for (int icalo = 0; icalo < n_RecCalo; icalo++) {

	int pindex = bRecCalo.getInt("pindex", icalo);
	int layer = bRecCalo.getInt("layer", icalo);

	// PCal, ECin, and ECout are part of the REC::Calorimeter Bank
	// they are distinguished by the variable "layer" 1 (PCal) 4 (ECin) and 7 (ECout)

	//  Part->ECindex relation
	if (layer == 1) {
	  indPCal[pindex] = icalo + 1;
	} else if (layer == 4) {
	  indECin[pindex] = icalo + 1;
	} else if (layer == 7) {
	  indECout[pindex] = icalo + 1;
	}

      }

      
      event.getStructure(bFowTgg);
      int nFT=bFowTgg.getRows();
      for(int iFT=0;iFT<nFT;iFT++){
	int pindex = bFowTgg.getInt("pindex"  ,iFT);
	int det    = bFowTgg.getInt("detector",iFT);
	if(det==10){
	  indFT[pindex]=iFT+1;
	}
      }
	 

      int CC_pidx;       
      event.getStructure(bRecCC);
      int nCC = bRecCC.getRows();
      for (int iCC = 0; iCC < nCC; iCC++) {
	int pindex = bRecCC.getInt("pindex", iCC);
	CC_pidx=pindex;
	int det = bRecCC.getInt("detector", iCC);
	// == Make sure the detector is HTCC
	if (det == DET_HTCC) {
	  // make the relation between Part->HTCCindex
	  indHTCC[pindex] = iCC + 1;
	}


      }


      int nSC = bRecSC.getRows();
      event.getStructure(bRecSC);
      for( int iSC = 0; iSC<nSC; iSC++){
	int pindex = bRecSC.getInt("pindex", iSC);
	int det    = bRecSC.getInt("detector",iSC);
	if(det==3 || det == 4) indSC[pindex] = iSC +1 ;
	//if(det==4){ indCTOF[pindex] = iSC+1;}
	//if(det==3){ indCND[pindex]  = iSC+1;}
      }
      

      event.getStructure(PART);



      
      event.getStructure(bBmtCr);
      int Bmt_ID=0,Bmt_Sector=0,Bmt_Layer=0,Bmt_Region=0;
      double nBmt = bBmtCr.getRows();
      //if(nBmt>0)bBmtCr.show();
      for (int iBmt  = 0; iBmt < nBmt; iBmt++) {

	Bmt_ID      = bBmtCr.getInt("ID",iBmt);
	Bmt_Sector  = bBmtCr.getInt("sector",iBmt);
	Bmt_Layer   = bBmtCr.getInt("layer",iBmt);
	Bmt_Region  = bBmtCr.getInt("region",iBmt);
  	
	
	h_BMT_ID       -> Fill(Bmt_ID);
	h_BMT_Region   -> Fill(Bmt_Region);
	h_BMT_sector   -> Fill(Bmt_Sector);
	
	
	
      }


      int CVT_ID=0,q_cvt=0,id_cvt=0,Cross_ID1=0,Cross_ID2=0,Cross_ID3=0,Cross_ID4=0,Cross_ID5=0,Cross_ID6=0,Cross_ID7=0,Cross_ID8=0,Cross_ID9=0;
      int tcvt=0;
      event.getStructure(bcvtTr);
      double ncvt = bcvtTr.getRows();
      for (int icvt  = 0; icvt < ncvt; icvt++) {
	CVT_ID     = bcvtTr.getInt("ID",icvt);
	Cross_ID1 = bcvtTr.getInt("Cross1_ID",icvt);
	Cross_ID2 = bcvtTr.getInt("Cross2_ID",icvt);
	Cross_ID3 = bcvtTr.getInt("Cross3_ID",icvt);
	Cross_ID4 = bcvtTr.getInt("Cross4_ID",icvt);
	Cross_ID5 = bcvtTr.getInt("Cross5_ID",icvt);
	Cross_ID6 = bcvtTr.getInt("Cross6_ID",icvt);
	Cross_ID7 = bcvtTr.getInt("Cross7_ID",icvt);
	Cross_ID8 = bcvtTr.getInt("Cross8_ID",icvt);
	Cross_ID9 = bcvtTr.getInt("Cross9_ID",icvt);
	int q_cvt    = bcvtTr.getInt("q",icvt);
	int id_cvt   = bcvtTr.getInt("ID",icvt);
	int FitMed   = bcvtTr.getInt("fittingMethod",icvt);
	double p_cvt = bcvtTr.getFloat("p",icvt);
	int pid_cvt  = bcvtTr.getInt("pid",icvt);

	indCVT[icvt]=icvt+1;
	CVT_FM1=FitMed;
	h_FitMed-> Fill(FitMed);
	

      }

      //if(npart>=1)PART.show();
     
      /////////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////
      double d1=0,d2=0,d3=0,d4=0;

      for (int ipart = 0; ipart < npart; ipart++) {
	
	int charge = PART.getInt("charge", ipart);  // Charge of the particle
	int stat   = PART.getInt("status", ipart);    // status of the particle
                    
	int  pid = PART.getInt  ("pid", ipart);      // PID of the particle
	float px = PART.getFloat("px" , ipart);      // X component of particl'es momentum 
	float py = PART.getFloat("py" , ipart);      // Y component of particl'es momentum 
	float pz = PART.getFloat("pz" , ipart);      // Z component of particl'es momentum
	float vz = PART.getFloat("vz" , ipart);      // Z component of particl'es vertex (cm)
	float vt = PART.getFloat("vt" , ipart);      // corrected vertex time 
	
	
	float chi2 = PART.getFloat("chi2pid", ipart); // Chi2 of assigned PID
	float beta = PART.getFloat("beta", ipart);   // particle beta  

	float  p   = sqrt(px * px + py * py + pz * pz);
	double th  = acos(pz / p) * TMath::RadToDeg();
	double phi = 30 + atan2(py/p,px/p)*TMath::RadToDeg() ;

	

	if (phi <   0) phi=phi+360;
	if (phi > 360) phi=phi-360;
	//	int sector = int(phi/60.) ;
	//	sector = sector + 1; 
	//h_phi -> Fill(phi);
	
	int  sector;
	double du_Pcal=0, dv_Pcal=0, dw_Pcal=0, du_ECin =0, dv_ECin =0, dw_ECin=0,du_ECout =0, dv_ECout =0, dw_ECout=0;
	double m2u_Pcal=0, m2v_Pcal=0, m2w_Pcal=0, m3u_Pcal =0, m3v_Pcal =0, m3w_Pcal=0;
	double m2u_ECin=0, m2v_ECin=0, m2w_ECin=0, m3u_ECin =0, m3v_ECin =0, m3w_ECin=0;
	double m2u_ECout=0,m2v_ECout=0,m2w_ECout=0,m3u_ECout=0, m3v_ECout=0, m3w_ECout=0;
	double lu_Pcal=0,lv_Pcal=0,lw_Pcal=0;
	double EPcal=0,EECin=0,EECout=0;
        float ECx=0,ECy=0;
	float time_Pcal,time_CC;
	double E_FT,x_FT,y_FT,z_FT,t_FT,stat_FT,r_FT,size_FT;
	double hx=0,hy=0,hz=0,x=0,y=0,z=0;

	


	int Cross_ID1=0,Cross_ID2=0,Cross_ID3=0,Cross_ID4=0,Cross_ID5=0,Cross_ID6=0,Cross_ID7=0,Cross_ID8=0,Cross_ID9=0;
	if(indCVT[ipart]>0){
	  Cross_ID1 = bcvtTr.getInt("Cross1_ID",indCVT[ipart] - Partconst);
	  Cross_ID2 = bcvtTr.getInt("Cross2_ID",indCVT[ipart] - Partconst);
	  Cross_ID3 = bcvtTr.getInt("Cross3_ID",indCVT[ipart] - Partconst);
	  Cross_ID4 = bcvtTr.getInt("Cross4_ID",indCVT[ipart] - Partconst);
	  Cross_ID5 = bcvtTr.getInt("Cross5_ID",indCVT[ipart] - Partconst);
	  Cross_ID6 = bcvtTr.getInt("Cross6_ID",indCVT[ipart] - Partconst);
	  Cross_ID7 = bcvtTr.getInt("Cross7_ID",indCVT[ipart] - Partconst);
	  Cross_ID8 = bcvtTr.getInt("Cross8_ID",indCVT[ipart] - Partconst);
	  Cross_ID9 = bcvtTr.getInt("Cross9_ID",indCVT[ipart] - Partconst);

	 
	
	
	   
	}
       

	// Getting energy depositions in PCal, ECin and ECout of the particle
	if (indPCal[ipart] > 0) {
	  EPcal      = bRecCalo.getFloat("energy", indPCal[ipart] - Partconst);
	  lu_Pcal    = bRecCalo.getFloat("lu"    , indPCal[ipart] - Partconst);
	  lv_Pcal    = bRecCalo.getFloat("lv"    , indPCal[ipart] - Partconst);
	  lw_Pcal    = bRecCalo.getFloat("lw"    , indPCal[ipart] - Partconst);
	  ECx        = bRecCalo.getFloat("x",      indPCal[ipart] - Partconst);
	  ECy        = bRecCalo.getFloat("y",      indPCal[ipart] - Partconst);
	  du_Pcal    = bRecCalo.getFloat("du",     indPCal[ipart] - Partconst);
	  dv_Pcal    = bRecCalo.getFloat("dv",     indPCal[ipart] - Partconst);
	  dw_Pcal    = bRecCalo.getFloat("dw",     indPCal[ipart] - Partconst);
	  m2u_Pcal   = bRecCalo.getFloat("m2u",    indPCal[ipart] - Partconst);  
	  m2v_Pcal   = bRecCalo.getFloat("m2v",    indPCal[ipart] - Partconst);
	  m2w_Pcal   = bRecCalo.getFloat("m2w",    indPCal[ipart] - Partconst);
	  m3u_Pcal   = bRecCalo.getFloat("m3u",    indPCal[ipart] - Partconst);
	  m3v_Pcal   = bRecCalo.getFloat("m3v",    indPCal[ipart] - Partconst);
	  m3w_Pcal   = bRecCalo.getFloat("m3w",    indPCal[ipart] - Partconst);
	  sector     = bRecCalo.getInt  ("sector", indPCal[ipart] - Partconst);
	  time_Pcal  = bRecCalo.getFloat("time",   indPCal[ipart] - Partconst);
	}
	if (indECin[ipart] > 0) {
	  EECin      = bRecCalo.getFloat("energy", indECin[ipart] - Partconst);
	  du_ECin    = bRecCalo.getFloat("du",     indECin[ipart] - Partconst);
	  dv_ECin    = bRecCalo.getFloat("dv",     indECin[ipart] - Partconst);
	  dw_ECin    = bRecCalo.getFloat("dw",     indECin[ipart] - Partconst);
	  m2u_ECin   = bRecCalo.getFloat("m2u",    indECin[ipart] - Partconst);
	  m2v_ECin   = bRecCalo.getFloat("m2v",    indECin[ipart] - Partconst);
	  m2w_ECin   = bRecCalo.getFloat("m2w",    indECin[ipart] - Partconst);
	  m3u_ECin   = bRecCalo.getFloat("m3u",    indECin[ipart] - Partconst);
	  m3v_ECin   = bRecCalo.getFloat("m3v",    indECin[ipart] - Partconst);
	  m3w_ECin   = bRecCalo.getFloat("m3w",    indECin[ipart] - Partconst);
	}
	if (indECout[ipart] > 0) {
	  EECout      = bRecCalo.getFloat("energy", indECout[ipart] - Partconst);
	  du_ECout    = bRecCalo.getFloat("du",     indECout[ipart] - Partconst);
	  dv_ECout    = bRecCalo.getFloat("dv",     indECout[ipart] - Partconst);
	  dw_ECout    = bRecCalo.getFloat("dw",     indECout[ipart] - Partconst);	  
	  m2u_ECout   = bRecCalo.getFloat("m2u",    indECout[ipart] - Partconst);
	  m2v_ECout   = bRecCalo.getFloat("m2v",    indECout[ipart] - Partconst);
	  m2w_ECout   = bRecCalo.getFloat("m2w",    indECout[ipart] - Partconst);
	  m3u_ECout   = bRecCalo.getFloat("m3u",    indECout[ipart] - Partconst);
	  m3v_ECout   = bRecCalo.getFloat("m3v",    indECout[ipart] - Partconst);
	  m3w_ECout   = bRecCalo.getFloat("m3w",    indECout[ipart] - Partconst);
	}

	double Etot = EPcal + EECin + EECout;
	double EEC  = EECin + EECout;

	float nphe = 0;
	 
	if (indHTCC[ipart] > 0) {
	  nphe = bRecCC.getFloat("nphe", indHTCC[ipart] - Partconst);
	  time_CC = bRecCC.getFloat("time", indHTCC[ipart] - Partconst);
	}

	if(indFT[ipart]>0){
	  E_FT = bFowTgg.getFloat("energy",    indFT[ipart] - Partconst);
	  x_FT = bFowTgg.getFloat("x",         indFT[ipart] - Partconst);
	  y_FT = bFowTgg.getFloat("y",         indFT[ipart] - Partconst);
	  z_FT = bFowTgg.getFloat("z",         indFT[ipart] - Partconst);
	  r_FT = bFowTgg.getFloat("radius",    indFT[ipart] - Partconst);
	  t_FT = bFowTgg.getFloat("time",      indFT[ipart] - Partconst);
	  size_FT = bFowTgg.getFloat("size",   indFT[ipart] - Partconst); 
	  stat_FT=bFowTgg.getFloat("status",   indFT[ipart] - Partconst);
	}

	if (indSC[ipart] > 0) {
	  hx = bRecSC.getFloat("hx", indSC[ipart] - Partconst);
	  hy = bRecSC.getFloat("hy", indSC[ipart] - Partconst);
	  hz = bRecSC.getFloat("hz", indSC[ipart] - Partconst);
	  x  = bRecSC.getFloat("x",  indSC[ipart] - Partconst);
	  y  = bRecSC.getFloat("y",  indSC[ipart] - Partconst);
	  z  = bRecSC.getFloat("z",  indSC[ipart] - Partconst);
	}
	
	//
	//
       	//
	//

	if (pid==11 && p>=2 && TMath::Abs(stat)>=2000 && TMath::Abs(stat)<4000  ) {
	  Electron_number++;  	  
	  Electron_sector[Electron_number]= sector;
	  P_el [Electron_number] = p;
	  P_elx[Electron_number] = px;
	  P_ely[Electron_number] = py;
	  P_elz[Electron_number] = pz;
	  vt_el[Electron_number] = vt;
	  el_phi  [Electron_number] = phi ;
	  el_theta[Electron_number] = th ;
	  //if(npart>=1)PART.show();
	}

	if (pid==-11 && p>=1.5 && TMath::Abs(stat)>=2000 && TMath::Abs(stat)<4000  ) {
	  Positron_number++;  	  
	  Positron_sector[Positron_number]= sector;
	  P_ps [Positron_number] = p;
	  P_psx[Positron_number] = px;
	  P_psy[Positron_number] = py;
	  P_psz[Positron_number] = pz;
	  vt_ps[Positron_number] = vt;
	  ps_phi  [Positron_number] = phi ;
	  ps_theta[Positron_number] = th ;
	}

	if(charge>0 && p>=0.4 &&  TMath::Abs(stat)>=2000 && TMath::Abs(stat)<4000 && TMath::Abs(chi2)<3 ){
	  Positive_FD_number++;
	  Positive_FD_sector[Positive_FD_number]= sector;
          P_FD_pos[Positive_FD_number]    = p;
	  P_FD_posx[Positive_FD_number]   = px;
	  P_FD_posy[Positive_FD_number]   = py;
	  P_FD_posz[Positive_FD_number]   = pz;
	  vt_FD_pos[Positive_FD_number]   = vt;
	  pos_FD_phi[Positive_FD_number]  = phi;
	  pos_FD_theta[Positive_FD_number]= th;

	  d1 = sqrt((x-hx)*(x-hx) + (y-hy)*(y-hy) + (z-hz)*(z-hz));
	  
	}
	
	if(charge<0 && pid!=11 && p>=0.4 && TMath::Abs(stat)>=2000 && TMath::Abs(stat)<4000 && TMath::Abs(chi2)<3){
	  Negative_FD_number++;
	  Negative_FD_sector[Negative_FD_number]= sector;
          P_FD_neg[Negative_FD_number]     = p;
	  P_FD_negx[Negative_FD_number]    = px;
	  P_FD_negy[Negative_FD_number]    = py;
	  P_FD_negz[Negative_FD_number]    = pz;
	  vt_FD_neg[Negative_FD_number]    = vt;
	  neg_FD_phi[Negative_FD_number]   = phi;
	  neg_FD_theta[Negative_FD_number] = th;
	  d2 = sqrt((x-hx)*(x-hx) + (y-hy)*(y-hy) + (z-hz)*(z-hz));
	}

	if(charge>0 /*&& pid==211*/ && p>=0.4 && TMath::Abs(stat)>=4000 && TMath::Abs(chi2)<1 ){
	  Positive_CD_number++;
          P_CD_pos[Positive_CD_number]     = p;
	  P_CD_posx[Positive_CD_number]    = px;
	  P_CD_posy[Positive_CD_number]    = py;
	  P_CD_posz[Positive_CD_number]    = pz;
	  vt_CD_pos[Positive_CD_number]    = vt;
	  pos_CD_phi[Positive_CD_number]   = phi;
	  pos_CD_theta[Positive_CD_number] = th;
	  d3 = sqrt((x-hx)*(x-hx) + (y-hy)*(y-hy) + (z-hz)*(z-hz));
	}
	
	if(charge<0 /*&& pid==-211*/ && pid!=11 && p>=0.4 && TMath::Abs(stat)>=4000 && TMath::Abs(chi2)<1 ){
	  Negative_CD_number++;
          P_CD_neg[Negative_CD_number]     = p;
	  P_CD_negx[Negative_CD_number]    = px;
	  P_CD_negy[Negative_CD_number]    = py;
	  P_CD_negz[Negative_CD_number]    = pz;
	  vt_CD_neg[Negative_CD_number]    = vt;
	  neg_CD_phi[Negative_CD_number]   = phi;
	  neg_CD_theta[Negative_CD_number] = th;
	  d4 = sqrt((x-hx)*(x-hx) + (y-hy)*(y-hy) + (z-hz)*(z-hz));
	}

	if(pid == -211 && p>=0.4 && TMath::Abs(stat)>=4000  ){ h_chi2_neg -> Fill(chi2);}
	if(pid ==  211 && p>=0.4 && TMath::Abs(stat)>=4000  ){ h_chi2_pos -> Fill(chi2);}
	

	 
      } // -----end of part loop-------

      //////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////
      
      if(Electron_number == 1 ) { el_counter++;
	if(Bmt_Sector==1){el_counter_s1++;
	  h_Electron_number_s1 -> Fill(Electron_number);}
	if(Bmt_Sector==2){el_counter_s2++;
	  h_Electron_number_s2 -> Fill(Electron_number);}
	if(Bmt_Sector==3){el_counter_s3++;
	  h_Electron_number_s3 -> Fill(Electron_number);}
	h_Electron_th     -> Fill(el_theta[1]);
	       
	if(Positive_CD_number>=1){pos_counter++;
	  h_Positive_CD_number -> Fill(Positive_CD_number);
	  for(int l=1;l<=Positive_CD_number;l++){
	    h_positive_angles->Fill(pos_CD_theta[l],pos_CD_phi[l]);
	    h_Positive_CD_th     -> Fill(pos_CD_theta[l]);}
	}
	if(Negative_CD_number>=1){neg_counter++;
	  h_Negative_CD_number -> Fill(Negative_CD_number);
	  for(int l=1;l<=Negative_CD_number;l++){h_negative_angles->Fill(neg_CD_theta[l],neg_CD_phi[l]);
	    h_Negative_CD_th     -> Fill(neg_CD_theta[l]);}
	}
	
	if(Electron_number == 1 && Positive_CD_number>=1 && Positive_FD_number>=1 && Electron_sector[1]!=Positive_FD_sector[Positive_FD_number]){
	  double d_phi = pos_CD_phi[1] - pos_FD_phi[1];
	  double d_th  = pos_CD_theta[1] - pos_FD_theta[1] ;
	  h_positive_angle_dp -> Fill(d_phi,d_th);
	  if(Bmt_Sector==1) h_positive_angle_dp_s1 -> Fill(d_phi,d_th);
	  if(Bmt_Sector==2) h_positive_angle_dp_s2 -> Fill(d_phi,d_th);
	  if(Bmt_Sector==3) h_positive_angle_dp_s3 -> Fill(d_phi,d_th);
		 
	  if(fabs(d_phi)<10 && fabs(d_th)<5){
	    pos_counter_cut++;
	    h_positive_angle_dp_cut_s1 -> Fill(d_phi,d_th);
	    if(Bmt_Sector==1){pos_counter_cut_s1++;
	      h_positive_angle_dp_cut_s1 -> Fill(d_phi,d_th);}
	    if(Bmt_Sector==2){pos_counter_cut_s2++;
	      h_positive_angle_dp_cut_s2 -> Fill(d_phi,d_th);}
	    if(Bmt_Sector==3){pos_counter_cut_s3++;
	      h_positive_angle_dp_cut_s3 -> Fill(d_phi,d_th);}
	  }
	}
	
	if(Electron_number == 1 && Negative_CD_number>=1 && Negative_FD_number>=1 && Electron_sector[1]!=Negative_FD_sector[Negative_FD_number]){
	  double d_phi = pos_CD_phi[1] - pos_FD_phi[1];
	  double d_th  = pos_CD_theta[1] - pos_FD_theta[1] ;
	  h_negative_angle_dp -> Fill(d_phi,d_th);
	  if(Bmt_Sector==1) h_negative_angle_dp_s1 -> Fill(d_phi,d_th);
	  if(Bmt_Sector==2) h_negative_angle_dp_s2 -> Fill(d_phi,d_th);
	  if(Bmt_Sector==3) h_negative_angle_dp_s3 -> Fill(d_phi,d_th);
	  
	  if(fabs(d_phi)<10 && fabs(d_th)<5){
	    neg_counter_cut++;
	    h_negative_angle_dp_cut -> Fill(d_phi,d_th);
	    if(Bmt_Sector==1){neg_counter_cut_s1++;
	      h_negative_angle_dp_cut_s1 -> Fill(d_phi,d_th);}
	    if(Bmt_Sector==2){neg_counter_cut_s2++;
	      h_negative_angle_dp_cut_s2 -> Fill(d_phi,d_th);}
	    if(Bmt_Sector==3){neg_counter_cut_s3++;
	      h_negative_angle_dp_cut_s3 -> Fill(d_phi,d_th);}
	  }
	}
      }
            
      ///////////////////////////////////
      ///////////////////////////////////


                
      //===============================================================================================================

      //      out_data<<test<<endl;  
      counter++;
    } // file event counter
  } // iFile counter
  
  
  out_data<<N_run<<'\t'<<el_counter<<'\t'<<pos_counter<<'\t'<<neg_counter<<'\t'<<pos_counter_cut<<'\t'<<neg_counter_cut<<'\t'<<50<<'\t'<<0.<<endl;
  
  gDirectory->Write();
  out_data.close();
  
  return 0;
}



