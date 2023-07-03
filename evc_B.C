#include <iostream>
#include <fstream>
#include <math.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <fstream>
#include <iostream>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLatex.h>

using namespace std;




void evc_B(){
 
   const double Lumin = 1/(3.3 * pow(10,39));
  TCanvas *c0 = new TCanvas("c0","c0",1280, 960);
  TStyle *style1 = new TStyle("Plain", "");
  style1->SetPalette(1); 
  gStyle->SetOptFit(1);
 
  
  style1->SetPalette(1);
  style1->SetOptFit();
  style1->cd();

  std::cout<<std::fixed;
  std::cout << std::setprecision(5);

  
  ifstream EF;
  ofstream test;
  EF.open("acc_MC.dat");
  test.open("test.dat");

  //double WW[]       = {4.03,4.09,4.15,4.21,4.27,4.33,4.39,4.45,4.51,4.57};
  // double WW_er[]    = {0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03};

  double NQ2[100],NW[100],Nt[100],Nacc[100],Nacce[100];
  int line=0;
  const int dot =8;
  double  Q2 ,W,t,E,bin,binvalue,MCbinvalue,RECbinvalue,errorbin;
   if(EF.is_open()) EF>>Q2>>W>>t>>E>>bin>>binvalue>>MCbinvalue>>RECbinvalue>>errorbin;
   
   while(1){line++;
     if(EF.eof()) break;

     NQ2[line] = Q2;
     NW[line] = W;
     Nt[line]= t ;
     //Nacc[line]=1/binvalue;
     //Nacce[line]=1/errorbin;
     Nacc[line]=MCbinvalue/RECbinvalue;
     Nacce[line]=(MCbinvalue/RECbinvalue)*sqrt(1/MCbinvalue + 1/RECbinvalue);
     

     test<<line<<'\t'<<NQ2[line]<<" "<<NW[line]<<" "<<Nt[line]<<" "<<Nacc[line]<<" "<<Nacce[line]<<endl;
     
     if(EF.is_open())EF>>Q2>>W>>t>>E>>bin>>binvalue>>MCbinvalue>>RECbinvalue>>errorbin;
   }
  int lmax = line; 
  double S[dot],Se[dot],SS[dot], NWW[dot],NWW_er[dot],Nk[dot]={0};

  int ni = 0 ,nk=0,nj=0;
  double step = 0.6/dot;
  
  for(int i=0;i<=lmax;i++){
    for( int k=0;k<dot;k++){
      if(NW[i]>=(4+step*k) && NW[i]<(4+step*(k+1))){Nk[k]++;
	S[k] = S[k]+Nacc[i];
	//Se[k]= Se[k]+Nacce[i];
	cout << "Count for S = " << k << ": " << S[k] << endl;
      }
    }
  }
  
  for( int j=0;j<dot;j++){
    SS[j] = Lumin*S[j];
    //cout<<j<<"  " <<S[j]<<endl;
    	NWW[j]= 4 + step*(j+0.5);
	NWW_er[j] = 0.5*step;
    //Se[j] =Lumin*Nk[j]; 
    //Se[j] = Lumin*(Nk[j]/S[j])*sqrt(1/Nk[j]+pow((S[j]/Se[j]),2));;
  }

  int ndot = dot -1 ; 
  auto *gr2 = new TGraph(ndot,NWW,SS);
  auto *gr3 = new TGraphErrors(ndot,NWW,SS,NWW_er,Se);
  //gr2->SetMarkerStyle(20);
  gr2->Draw("AP");
   gr2 -> Print();
   
   /*

  double Sigma[]    = {Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10};
  double Sigma_er[] = {Y1e,Y2e,Y3e,Y4e,Y5e,Y6e,Y7e,Y8e,Y9e,Y10e};
  double nul[]      = {0,0,0,0,0,0,0,0,0,0};
  double WW[]       = {4.03,4.09,4.15,4.21,4.27,4.33,4.39,4.45,4.51,4.57};
  double WW_er[]    = {0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03};

  //double WW[]       = {8.44,8.92,9.4,9.88,10.36};
  //double WW_er[]    = {0.24,0.24,0.24,0.24,0.24};


  
  auto *gr1 = new TGraphErrors(10,WW,Sigma,WW_er,Sigma_er);
  //auto *gr1 = new TGraph(4,WW,Sigma);

  //double maxY = gr1->GetYaxis()->GetXmax();
  //double newMaxY = maxY*1000;
  //gr1->GetYaxis()->SetRangeUser(0., newMaxY);
  
  cout<<Y1<<endl;
  cout<<Y2<<endl;
  cout<<Y3<<endl;
  cout<<Y4<<endl;
  cout<<Qtot<<endl;
  cout<<charge_e<<endl;
  cout<<Number_pr<<endl;
  cout<<Lumin<<endl;
  //cout<<i1<<" "<<i2<<" "<<i3<<" "<<i4<<" "<<i5<<endl;
  

  // gr1->SetLineWidth(2);
  // //gr1->SetLineStyle(9);
  // gr1->SetMarkerColor(2);
  // gr1->SetMarkerStyle(20);
  
      
  //gr1->SetTitle("positive cases sector = 3");
  // gr1->Fit("f_pol_1", "GeV","",  5 ,50.5);
  //gr1->GetXaxis()->SetTitle("W [GeV]");
  //gr1->GetYaxis()->SetTitle("Sigma");
  
  
  gr1->Draw("AP");
  gr1 -> Print();
  //cout<<"####################### "<<f_pol_1_0<<endl;
  // cout<<"####################### "<<f_pol_2_0<<endl;
  */ 
  test.close();

  gDirectory -> Write();
   
}

