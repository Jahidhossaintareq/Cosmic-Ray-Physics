/*****************************************************************************  
 *                          jahid_hossain_1.C                                * 
 *        Macro to analyze the TREE created with anacorsika_1.cc             *
 *        Version : 5/Jun/2013                                               *
 *                                                                           * 
 *        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        *
 *        %                This code is modified by                 %        * 
 *        %                   JAHID HOSSAIN                         %        *
 *        %         Department of Theoretical Physics               %        *
 *        %                 University of Dhaka                     %        *
 *        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        *
 *****************************************************************************/

#include <TH1>
#include <TFile>
#include <TLegend>
#include <TCanvas>
#include <TMath>

void jahid_hossain_1() {

gROOT->Reset();
gStyle->SetOptStat("nemr");

Char_t infilnam[200];
Char_t outfil[200]; 
Char_t outplot[200];  
FILE *wrtout;

const Int_t maxmu=5000;  // maximum number of muons forseen in one event

Int_t ias;
Int_t num_ev_run, nmuev,nmuevx,nmuevy,nmuinali;
Float_t eneprim; 
Int_t muplus[maxmu];
Int_t muminus[maxmu];
Float_t emuev[maxmu];
Float_t tetamuev[maxmu];
Float_t phimuev[maxmu];

// %%%%%%%% Number of file to be analysed and output connected (sufix) %%%%%%%%%


Int_t sufix=1;
sprintf(outfil,"ptest_ali_31015_1016_%d.out",sufix);
sprintf(outplot,"ptest_ali_31015_1016_%d.root",sufix);


// %%%%%%%%%%%%%%%%%%%%%%%%  Open file to write infos  %%%%%%%%%%%%%%%%%%%%%%%%%


wrtout = fopen(outfil,"wt");


// %%%%%%%%%%%%%%%%%%% Open root file in which save the plots %%%%%%%%%%%%%%%%%%


//TFile *saveplot = new TFile("p7350_1016_1018_1year_1.root","RECREATE","plot_tree");
TFile *saveplot = new TFile(outplot,"RECREATE","plot_tree");


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



// *************************   Energy of the Primary ***************************

TH1F *heneprim=new TH1F("Eneprim #Epsilon [PeV] ","Proton Energy Spectrum for E = 0.1 PeV-3.0 PeV ;Primary Energy #Epsilon [PeV]; Num.of events",100,0.0,3.0);

// ************************ Muon Multiplicity Distribution ********************

 TH1F *h2=new TH1F("Around ALICE","Mult. Dist. around ALICE, Nmu>0 (#Epsilon_{p}= 0.1PeV - 3.0PeV); Number of muons #Nu_{#mu}; Number of events",20,0,20);
//TH1F *h6=new TH1F("Inside ALICE","Muon Mult. Dist. inside ALICE Nmu>0 ; Number of muons; Number of events",20,0.0,20);

// *************************  Muon Energy Distributions ************************

TH1F *hemuev=new TH1F("#Nu_{#mu} vs #Epsilon_{#mu} [GeV]","Muon Energy Spectrum (#Epsilon_{p}=0.1PeV - 3.0PeV);Muon Energy #Epsilon_{#mu} [GeV];Number of Muon #Nu_{#mu}",100,0.,1000.);

TH1F *hemuevp=new TH1F("R vs #Epsilon_{#mu} ","Muon Energy Distribution;Muon Energy #Epsilon_{#mu} [GeV];#mu^{+}",10,0.,150.);

TH1F *hemuevm=new TH1F("Mu Energy Dist.[GeV]","Mu Energy Distribution;Muon Energy #Epsilon_{#mu} [GeV];#mu^{-}",10,0.,150.);
 
TH1F *hecostatp=new TH1F("R vs #Epsilon_{#mu} cos#theta [GeV]"," ;#Epsilon_{#mu} cos#theta [GeV];#mu^{+}",15,0.,150.);

TH1F *hecostatm=new TH1F("R vs #Epsilon_{#mu} cos#theta"," ;#Epsilon_{#mu} cos#theta [GeV];#mu^{-}",15,0.,150.);

// ************************  Zenithal Distribution ******************************

TH1F *htetanmu=new TH1F("#Nu_{#mu} vs #theta","Zenithal angle distribution;Zenith angle #theta[deg] ; Number of Muon #Nu_{#mu}",100,0.,50.);

TH1F *htetanmup=new TH1F("R vs #theta ","Theta Primary vs Nmu;Zenith angle #theta[deg] ; Number of #mu+",10,0.,50.);

TH1F *htetanmum=new TH1F("Teta primary vs Nmu","Theta Primary vs Nmu;Zenith angle #theta[deg] ; Number of #mu-",10,0.,50.);

// ******************* Azimuthal Distribution **********************************

TH1F *hphinmu=new TH1F("#Nu_{#mu} vs #phi","Azimuthal angle distribution;Azimuthal angle #phi[deg] ; Number of Muon #Nu_{#mu}",20,0.,360.);

TH1F *hphinmup=new TH1F("R vs #phi ","Phi Primary vs Nmu;Azimuthal angle #phi[deg] ; Number of #mu+",10,0.,360.);

TH1F *hphinmum=new TH1F("Phi primary vs Nmu","Phi Primary vs Nmu;Azimuthal angle #phi[deg] ; Number of #mu-",10,0.,360.);

// ****************************** Muon plus ************************************

TH1F *hmuplus=new TH1F("#Nu_{#mu^{+}}","Muon Charge Distribution (#Epsilon_{p}= 0.1PeV - 3.0PeV);Charged Muon #mu^{+} & #mu^{-};Number of Muons #Nu_{#mu}",10,-5.,5.); 

// ***************************** Muon minus ************************************

TH1F *hmuminus=new TH1F("#Nu_{#mu^{-}}","Muon Charge Distribution;Charged Muon #mu^{+} & #mu^{-}; #Nu_{#mu}",10,-5.,5.); 

// ************************** Mu Charge Ratio caclulation ***********************

TH1F *hratiot=new TH1F("R vs #theta","Muon Charge Ratio Distribution (#Epsilon_{p}= 0.1PeV - 3.0PeV);Zenith angle #theta[deg];Charge Ratio R = #mu^{+}/#mu^{-}",10,0.,50.);

TH1F *hratiop=new TH1F("R vs #phi","Muon Charge Ratio Distribution (#Epsilon_{p}= 0.1PeV - 3.0PeV);Azimuthal angle #phi[deg];Charge Ratio R = #mu^{+}/#mu^{-}",10,0.,360.);

TH1F *hratio=new TH1F("R vs #Epsilon_{#mu}","Muon Charge Ratio Distribution (#Epsilon_{p}= 0.1PeV - 3.0PeV);Muon Energy #Epsilon_{#mu} [GeV];Charge Ratio R = #mu^{+}/#mu^{-}",10,0.,150.);
TH1F *hratioc=new TH1F("R vs #Epsilon_{#mu}cos#theta","Muon Charge Ratio Distribution (#Epsilon_{p}= 0.1PeV - 3.0PeV);Muon Energy  #Epsilon_{#mu}cos#theta  [GeV];Charge Ratio R = #mu^{+}/#mu^{-}",15,0.,150.);

// %%%%%%%%%%%%%%%%%%%%%%%%%%% Open the root files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


const Int_t nfile1=1;   // first root 
const Int_t nfile2=1;   // last root 
Int_t iloopf;
Int_t feentries;
Int_t iloopf2;

for(iloopf=nfile1;iloopf<nfile2+1;iloopf++){

if(iloopf==1){
sprintf(infilnam,"ptest_tree_31015_1016_%d.root",iloopf);
TFile *ffe = new TFile(infilnam);
TTree *tfe = (TTree*)ffe->Get("AliCorsEv");
Int_t maxnumev=10000000;  // maximum number of event to be used
feentries = (Int_t)tfe->GetEntries();
if(maxnumev<feentries)feentries=maxnumev;  // max.num.of events for this root file
fprintf(wrtout,"%d) === OPEN FILE : %s \n",iloopf,infilnam);
fprintf(wrtout," Entries =  %d    Max.Num.Ev. = %d \n \n",feentries,maxnumev);
cout << iloopf << ")  OPEN FILE : " << infilnam << "  Entries = " << feentries << "   Max.num.ev. = " <<  maxnumev << endl;
}

if(iloopf==2){
sprintf(infilnam,"ptest_tree_1016_31016_%d.root",iloopf);
TFile *ffe = new TFile(infilnam);
TTree *tfe = (TTree*)ffe->Get("AliCorsEv");
Int_t maxnumev=10000000;  // maximum number of event to be used
feentries = (Int_t)tfe->GetEntries();
if(maxnumev<feentries)feentries=maxnumev;  // max.num.of events for this root file
fprintf(wrtout,"%d) === OPEN FILE : %s \n",iloopf,infilnam);
fprintf(wrtout," Entries =  %d    Max.Num.Ev. = %d \n \n",feentries,maxnumev);
cout << iloopf << ")  OPEN FILE : " << infilnam << "  Entries = " << feentries << "   Max.num.ev. = " <<  maxnumev << endl;
}


// %%%%%%%%%%%%%%%%%%%%%%%%% SET TREE BRANCH ADDRESS %%%%%%%%%%%%%%%%%%%%%%%%%%%

tfe->SetBranchAddress("eneprim",&eneprim);
//tfe->SetBranchAddress("nmuinali",&nmuinali);
tfe->SetBranchAddress("nmuev",&nmuev);
tfe->SetBranchAddress("emuev",&emuev);
tfe->SetBranchAddress("tetamuev",&tetamuev);
tfe->SetBranchAddress("phimuev",&phimuev);
tfe->SetBranchAddress("nmuevx",&nmuevx);
tfe->SetBranchAddress("nmuevy",&nmuevy);
tfe->SetBranchAddress("muplus",&muplus);
tfe->SetBranchAddress("muminus",&muminus);

// %%%%%%%%%%%%%%%%%%%%%%%%% START LOOP OVER EVENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%


for(Int_t feentry=0;feentry<feentries;feentry++){

      tfe->GetEntry(feentry);
 
      heneprim->Fill(eneprim/1000000); 

       if(nmuev>0){
  	  h2->Fill(nmuev);  // muon mult. distr. Nmu>0 around ALICE
		} 

       //      if(nmuinali>0){
       // h6->Fill(nmuinali);  // muon mult. distr. Nmu>0 inside Alice
       //	} 

for(ias=0;ias<nmuev;ias++){
                    
      hemuev->Fill(emuev[ias]);                                             
      htetanmu->Fill(tetamuev[ias]);
      hphinmu->Fill(phimuev[ias]);
}
  
for(ias=0;ias<nmuevx;ias++){
                    
      hemuevp->Fill(emuev[ias]);
      hmuplus->Fill(muplus[ias]);                                             
      htetanmup->Fill(tetamuev[ias]);  
      hphinmup->Fill(phimuev[ias]);
      hecostatp->Fill(emuev[ias] * cos(tetamuev[ias]));
  } 
    
for(ias=0;ias<nmuevy;ias++){
      
       hemuevm->Fill(emuev[ias]);
       hmuminus->Fill(muminus[ias]);
       htetanmum->Fill(tetamuev[ias]);
       hphinmum->Fill(phimuev[ias]);       
       hecostatm->Fill(emuev[ias]*cos(tetamuev[ias]));
 }

} 

// %%%%%%%%%%%%%%%%%%%%%%%%%%% End of feentry loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

  cout << " Terminated loop over events -->  feentry = " << feentry << endl << endl;
//    Close the input file
  ffe->Close();
 }  


// %%%%%%%%%%%%%%%%%%%%%%%%%%% END LOOP ON ROOT FILES %%%%%%%%%%%%%%%%%%%%%%%%%%


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START  PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


// ******************* Plot of Primary Energy Spectrum *************************

TCanvas *a1=new TCanvas("a1","Primary Energy Spectrum",200,10,700,500);
  a1->cd(1);
  gStyle->SetOptStat(1000000211);
  heneprim->SetLineColor(kBlue+1);
  heneprim->SetLineWidth(2);
  heneprim->SetFillColor(kSpring-1);
  heneprim->GetYaxis()->SetTitleSize(20);
  heneprim->GetYaxis()->SetTitleFont(43);
  heneprim->GetYaxis()->SetTitleColor(kRed);
  heneprim->GetYaxis()->SetLabelColor(kRed); 
  heneprim->GetXaxis()->SetTitleSize(20);
  heneprim->GetXaxis()->SetTitleFont(43);
  heneprim->GetXaxis()->SetTitleColor(kRed);
  heneprim->GetXaxis()->SetLabelColor(kRed);    
  heneprim->Draw(); 

// *********************** plot of Muon multiplicity distribution **************

  TCanvas *m1=new TCanvas("m1","Muons Multiplicity Distribution",200,10,700,500);

  //m1->Divide(2,1);
  m1->cd(1);
  h2->SetLineColor(kBlack);
  h2->SetLineWidth(2);
  h2->SetFillColor(kBlue+1);
  h2->GetYaxis()->SetTitleSize(20);
  h2->GetYaxis()->SetTitleFont(43);
  h2->GetYaxis()->SetTitleColor(kRed);
  h2->GetYaxis()->SetLabelColor(kRed); 
  h2->GetXaxis()->SetTitleSize(20);
  h2->GetXaxis()->SetTitleFont(43);
  h2->GetXaxis()->SetTitleColor(kRed);
  h2->GetXaxis()->SetLabelColor(kRed);    
  h2->Draw();
  /* m1->cd(2);
  h6->SetLineColor(kBlue+1);
  h6->SetLineWidth(2);
  h6->SetFillColor(kGray+2);
  h6->GetYaxis()->SetTitleSize(20);
  h6->GetYaxis()->SetTitleFont(43);
  h6->GetYaxis()->SetTitleColor(kRed);
  h6->GetYaxis()->SetLabelColor(kRed); 
  h6->GetXaxis()->SetTitleSize(20);
  h6->GetXaxis()->SetTitleFont(43);
  h6->GetXaxis()->SetTitleColor(kRed);
  h6->GetXaxis()->SetLabelColor(kRed);    
  h6->Draw(); 
  */

// ************************ Plot of Charged muon vs Nmu ************************

TCanvas *b1=new TCanvas("b1","Charged Muon",200,10,700,500);
  b1->cd(1);
  gStyle->SetOptStat(1000000011);
  hmuplus->SetLineColor(kBlack);
  hmuplus->SetLineWidth(2);  
  hmuplus->SetFillColor(kAzure-1);  
  hmuplus->GetYaxis()->SetTitleSize(20);
  hmuplus->GetYaxis()->SetTitleFont(43);
  hmuplus->GetYaxis()->SetTitleColor(kRed);
  hmuplus->GetYaxis()->SetLabelColor(kRed); 
  hmuplus->GetXaxis()->SetTitleSize(20);
  hmuplus->GetXaxis()->SetTitleFont(43);
  hmuplus->GetXaxis()->SetTitleColor(kRed);
  hmuplus->GetXaxis()->SetLabelColor(kRed); 
  hmuplus->Draw();
  b1->Update();
  gStyle->SetOptStat(1000000011);
  hmuminus->SetLineColor(kBlue);
  hmuminus->SetLineWidth(2);
  hmuminus->SetFillColor(kGray+2);  
  hmuminus->Draw("sames");

  TLegend *feleg = new TLegend(0.7,0.9,0.9,0.7);
  feleg->SetTextFont(70);
  feleg->SetTextSize(0.030);
  feleg->SetTextColor(kGray+3);
  feleg->AddEntry(hmuplus,"Mu+");
  feleg->AddEntry(hmuminus,"Mu-"); 
  feleg->Draw();
  
// *********************** Plot of Nmu vs muon energy **************************

TCanvas *c1=new TCanvas("c1","Muon Energy Spectrum",200,10,700,500);
 c1->cd(1);
 hemuev->SetLineColor(kBlack);
 hemuev->SetFillColor(kGray+1);
 hemuev->SetLineWidth(2);
 hemuev->GetYaxis()->SetTitleSize(20);
 hemuev->GetYaxis()->SetTitleFont(43);
 hemuev->GetYaxis()->SetTitleColor(kRed);
 hemuev->GetYaxis()->SetLabelColor(kRed); 
 hemuev->GetXaxis()->SetTitleSize(20);
 hemuev->GetXaxis()->SetTitleFont(43);
 hemuev->GetXaxis()->SetTitleColor(kRed);
 hemuev->GetXaxis()->SetLabelColor(kRed);
 //hemuev->Fit("landau"); 
 hemuev->Draw();

 /*   
// *********************** Plot of Nmu vs zenith angle **********************

TCanvas *d1=new TCanvas("d1","zenithal Angle Distribution",200,10,700,500);
 d1->Divide(1,2); 
 d1->cd(1);
 htetanmu->SetMarkerColor(kBlue+2); 
 htetanmu->SetMarkerStyle(22);
 htetanmu->SetMarkerSize(1.5);
 htetanmu->GetYaxis()->SetTitleSize(20);
 htetanmu->GetYaxis()->SetTitleFont(43);
 htetanmu->GetYaxis()->SetTitleColor(kRed);
 htetanmu->GetYaxis()->SetLabelColor(kRed); 
 htetanmu->GetXaxis()->SetTitleSize(20);
 htetanmu->GetXaxis()->SetTitleFont(43);
 htetanmu->GetXaxis()->SetTitleColor(kRed);
 htetanmu->GetXaxis()->SetLabelColor(kRed);
 htetanmu->Fit("pol4"); 
 htetanmu->Draw("E");
 d1->cd(2);
 hphinmu->SetMarkerColor(kRed); 
 hphinmu->SetMarkerStyle(22);
 hphinmu->SetMarkerSize(1.5);
 hphinmu->GetYaxis()->SetTitleSize(20);
 hphinmu->GetYaxis()->SetTitleFont(43);
 hphinmu->GetYaxis()->SetTitleColor(kRed);
 hphinmu->GetYaxis()->SetLabelColor(kRed); 
 hphinmu->GetXaxis()->SetTitleSize(20);
 hphinmu->GetXaxis()->SetTitleFont(43);
 hphinmu->GetXaxis()->SetTitleColor(kRed);
 hphinmu->GetXaxis()->SetLabelColor(kRed);
 // hphinmu->Fit("pol4"); 
 hphinmu->Draw("E");
 */

// *********************** Plot of ratio vs zenith angle ******************* 

TCanvas *r1=new TCanvas("r1"," Mu Charge Ratio vs zenith angle",200,10,1000,500);

  r1->Divide(2,1);
  r1->cd(1);

  hratiot->SetLineColor(kBlack);
  hratiot->SetLineWidth(2); 
  hratiot->SetMarkerColor(kBlue+2); 
  hratiot->SetMarkerStyle(22);
  hratiot->SetMarkerSize(1.5);   
  hratiot->SetStats(0);  
  hratiot->Sumw2();
  hratiot->Divide(htetanmup,htetanmum);
  hratiot->GetYaxis()->SetTitleSize(20);
  hratiot->GetYaxis()->SetTitleFont(43);
  hratiot->GetYaxis()->SetTitleColor(kRed);
  hratiot->GetYaxis()->SetLabelColor(kRed); 
  hratiot->GetXaxis()->SetTitleSize(20);
  hratiot->GetXaxis()->SetTitleFont(43);
  hratiot->GetXaxis()->SetTitleColor(kRed);
  hratiot->GetXaxis()->SetLabelColor(kRed); 
  // hratiot->Fit("pol3");
  hratiot->Draw("E1"); 

  // ******************** plot of charge ratio vs azimuthal angle *************


  //TCanvas *p1=new TCanvas("p1"," Mu Charge Ratio vs Azimuthal angle",200,10,700,500);
  // p1->cd(1);
  r1-> cd(2);
  hratiop->SetLineColor(kBlack);
  hratiop->SetLineWidth(2); 
  hratiop->SetMarkerColor(kRed); 
  hratiop->SetMarkerStyle(22);
  hratiop->SetMarkerSize(1.5);   
  hratiop->SetStats(0);  
  hratiop->Sumw2();
  hratiop->Divide(hphinmup,hphinmum);
  hratiop->GetYaxis()->SetTitleSize(20);
  hratiop->GetYaxis()->SetTitleFont(43);
  hratiop->GetYaxis()->SetTitleColor(kRed);
  hratiop->GetYaxis()->SetLabelColor(kRed); 
  hratiop->GetXaxis()->SetTitleSize(20);
  hratiop->GetXaxis()->SetTitleFont(43);
  hratiop->GetXaxis()->SetTitleColor(kRed);
  hratiop->GetXaxis()->SetLabelColor(kRed); 
  // hratiop->Fit("pol3");
  hratiop->Draw("E1"); 



// ********************* Plot of Mu+ vs Muon energy ****************************

TCanvas *e2=new TCanvas("e2","Mu+ vs Muon Energy",200,10,700,500);
  e2->cd(1);
 
  hemuevp->SetLineColor(kBlue+1);
  hemuevp->SetLineWidth(2);
  hemuevp->SetFillColor(kSpring-1); 
  hemuevp->Draw(); 
  hemuevp->Print("all");

// ******************** Plot of Mu- vs Muon Energy *****************************

TCanvas *e3=new TCanvas("e3","Mu- vs Muon Energy",200,10,700,500);
  e3->cd(1);
 
  hemuevm->SetLineColor(kBlue+1);
  hemuevm->SetLineWidth(2);
  hemuevm->SetFillColor(kRed); 
  hemuevm->Draw(); 
  hemuevm->Print("all");

// ****************** Plot of Muon Charge ratio vs Muon energy **************** 

TCanvas *r2=new TCanvas("r2","Mu Charge Ratio vs Mu Energy",200,10,700,500);
  r2->cd(1);

  hratio->SetLineColor(kBlack);
  hratio->SetLineWidth(2);
  hratio->SetMarkerColor(kBlue); 
  hratio->SetMarkerStyle(22);
  hratio->SetMarkerSize(1.5); 
  hratio->SetStats(0);  
  hratio->Sumw2();  
  hratio->Divide(hemuevp,hemuevm);
  hratio->GetYaxis()->SetTitleSize(20);
  hratio->GetYaxis()->SetTitleFont(43);
  hratio->GetYaxis()->SetTitleColor(kRed);
  hratio->GetYaxis()->SetLabelColor(kRed); 
  hratio->GetXaxis()->SetTitleSize(20);
  hratio->GetXaxis()->SetTitleFont(43);
  hratio->GetXaxis()->SetTitleColor(kRed);
  hratio->GetXaxis()->SetLabelColor(kRed); 
  // hratio->Fit("pol8");  
  hratio->Draw("E1"); 
  hratio->Print("all");


// **************** Plot of Charg ratio vs Energy*cos(theta) *******************

   
//TCanvas *r3=new TCanvas("r3","Charge Ratio vs E*cos(theta)",200,10,700,500);
  r3->cd(1);

  hratioc->SetLineColor(kBlack);
  hratioc->SetLineWidth(2);
  hratioc->SetMarkerColor(kBlue); 
  hratioc->SetMarkerStyle(22);
  hratioc->SetMarkerSize(1.5);
  hratioc->SetStats(0); 
  hratioc->Sumw2();
  hratioc->Divide(hecostatp,hecostatm);
  hratioc->GetYaxis()->SetTitleSize(20);
  hratioc->GetYaxis()->SetTitleFont(43);
  hratioc->GetYaxis()->SetTitleColor(kRed);
  hratioc->GetYaxis()->SetLabelColor(kRed); 
  hratioc->GetXaxis()->SetTitleSize(20);
  hratioc->GetXaxis()->SetTitleFont(43);
  hratioc->GetXaxis()->SetTitleColor(kRed);
  hratioc->GetXaxis()->SetLabelColor(kRed); 
  //  hratioc->Fit("pol8");  
  // hratioc->Rebin();
  hratioc->Draw("E1"); 
  

saveplot->Write();
//saveplot->Close();

}
 






















