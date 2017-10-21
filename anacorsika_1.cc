/*******************************************************************************
*                                                                              *
*      Cosmic Program :                                                       *
* Read  Corsika files, analyze and create root file as input                 * 
* for further analysis with Montecarlo                                      *    
*                                                                           *   
*     Version: 1/June/2015                                                   * 
*     anacorsika_1.cc                                                         *
*                                                                              *
* Sampling Area 205 x 205 m**2 above Alice, the core in taken                  *
* randomly in this area                                                       *
* Alice underground 40 m is a rectangle 3 x 5 m**2 located in the center     * 
* and the energy threshold is E > Ethr GeV/cos(teta) at surface             *  
* To prepare the input for AliRoot all the muons arriving in an             *  
* area 6 x 6 m**2 centered in Alice are written in the tree                  * 
* The E Threshold for these muons is 1 GeV at surface                         *
*                                                                              *
*   To compile use :                                                           *
*   make -f make_corsika  creates anacorsika_1.exe                            * 
*   To run  :                                                                * 
*   ./anacorsika_1.exe  >  outputfile                                       *  
*                                                                            * 
*         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         *
*         %                This code is modified by                 %          *
*         %                   JAHID HOSSAIN                         %          *
*         %         Department of Theoretical Physics               %         *
*         %                 University of Dhaka                     %        *
*         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        *
*******************************************************************************/

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include "aliana_corsika.h"
using namespace std;

extern void AnaEv(int n);    

// *****************************************************************************

Int_t nn;
Float_t run,slope,nev,part_id,eneprim,emin,emax;
Float_t eneprimPeV,enerefPeV;     //ak
Float_t gammagen,gammaknee,deltagenknee,deltakneegen;   //ak
Int_t tagg30;     //ak
Float_t rapene;    //ak
Int_t num_ev_tot,num_ev_run,num_ev_ene;
Int_t num_ev_tot_g30,num_ev_run_g30,num_ev_ene_g30;   //ak
Int_t num_ev_test;

Int_t nmutdis[LOOPCORE];
Float_t primev;
Float_t nmu[ILOOP],tnmu[ILOOP],tdensmu[ILOOP],emur[ILOOP],temur[ILOOP];
Float_t tnmu_app_zero[ILOOP],tnmu_core_zero[ILOOP];
Float_t tnmu_app_ecut[ILOOP],tnmu_core_ecut[ILOOP];
Float_t nmutot,sizenkg,nmut15;
Float_t rag[ILOOP+1],area[ILOOP];
Float_t media_nmutot,media_nmut15,media_sizenkg;
Float_t err_media_nmutot,err_media_nmut15;

Float_t ltot=20500;      // total lenght in X and Y of the surface above Alice in cm.
const Int_t maxmu=4000;  // max. number of muons forseen inside large area around Alice
Float_t nmuali,nmualix,nmualiy;          // number of muons inside large area around Alice (6x6 m**2)

Int_t idchmu[maxmu];
Int_t muplus[maxmu];
Int_t muminus[maxmu];

Int_t   tagmucell[maxmu];   
Float_t tetamucell[maxmu];
Float_t phimucell[maxmu];
Float_t emucell[maxmu];

Int_t   tagmu[maxmu];   
Float_t tetamuev[maxmu];
Float_t phimuev[maxmu];
Float_t emuev[maxmu];

Int_t idcell,nmuev,nmuvec,nmuevx,nmuvecx,nmuevy,nmuvecy;
Float_t xco,yco;                      // x core, y core 
Int_t nmuinali,nmuinalix,nmuinaliy;   // number of muons inside Alice (5 x 3 m**2)
 
Float_t a1;
Float_t rea[50];
Int_t i;

union heads{
  Float_t header ;
  Char_t headname[4] ;
}   head ;

Float_t firstlast ;
FILE *bufana;

Char_t infilnam[150];

// *****************************************************************************

Int_t main(){

  Int_t loopf; 
  Int_t loopf_ene;    //ak
  
//   Index gamma of the generation (CORSIKA) and our gamma of the knee
  gammagen = -2.7;   //ak  gamma used in Corsika
  gammaknee = -3.0;   //ak gamma used after the knee
  deltagenknee = gammagen-gammaknee;      //ak
  deltakneegen = gammaknee-gammagen;      //ak
  
  //  Input for analyze some simulated runs  
  const Int_t nfile1=1;  // first file to be analysed
  const Int_t nfile2=1;  // last file to be analysed 
  const Int_t nfilevec=30;   // maximum number of files
  Int_t nevf[nfilevec];  //max. number of events per file
  Int_t nevf_ene[5]; //ak max num. of events for each range of energy


    nevf_ene[0]=1000000; // ak  gamma=-3.0  3 10**15-10**16  files:0-5
    nevf_ene[1]=0;  // ak  10**16-3 10**16          files:6-9
    nevf_ene[2]=0;   // ak  3 10**16-10**17         files:10-13
    nevf_ene[3]=0;    // ak  10**17-3 10**17        files:14-15
    nevf_ene[4]=0;     // ak  3 10**17-10**18       files:16
  
  for(i=0; i<nfilevec; i++){
    nevf[i]=10000000;
  }
  

//  %%%%%%%%%%%%%%%%%%%%%%%%%%% Prepare the TREE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

  TFile *ftree1 = new TFile("ptest_tree_31015_1016_1.root","RECREATE","alicors_tree");  
  TTree *ev_tree = new TTree("AliCorsEv","AliCors_ev_tree");

  ev_tree->Branch("eneprim",&eneprim,"eneprim/F");  
  ev_tree->Branch("nmuev",&nmuev,"nmuev/I"); 
  ev_tree->Branch("nmuevx",&nmuevx,"nmuevx/I"); 
  ev_tree->Branch("nmuevy",&nmuevy,"nmuevy/I");
  ev_tree->Branch("nmuinali",&nmuinali,"nmuinali/I");   
  ev_tree->Branch("emuev",emuev,"emuev[nmuev]/F");
  ev_tree->Branch("tetamuev",tetamuev,"tetamuev[nmuev]/F");  
  ev_tree->Branch("phimuev",phimuev,"phimuev[nmuev]/F");
  ev_tree->Branch("idchmu",&idchmu,"idchmu[nmuev]/I");
  ev_tree->Branch("muplus",&muplus,"muplus[nmuevx]/I");
  ev_tree->Branch("muminus",&muminus,"muminus[nmuevy]/I");
 

// %%%%%%%%%%%%%%%%%%%%%%% Some initializations for counting %%%%%%%%%%%%%%%%%%%

  num_ev_tot = 0;
  num_ev_run = 0;
  num_ev_test = 0;
  num_ev_ene = 0;
  num_ev_tot_g30 = 0;  //ak
  num_ev_run_g30 = 0;  //ak
  num_ev_ene_g30 = 0;  //ak                                
  media_nmutot  = 0 ;
  media_nmut15  = 0 ;
  media_sizenkg = 0 ;

// *********** Initialization for lateral density of muons [in cm] *************

  for(i=0; i<=ILOOP; i++){
    rag[i]= 200*i;
  }
  for(i=0; i<ILOOP; i++){
    area[i]=PIGR * ((rag[i+1]*rag[i+1])-(rag[i]*rag[i]));
    tnmu_app_zero[i] = tnmu_core_zero[i] = 0;
    tnmu_app_ecut[i] = tnmu_core_ecut[i] = 0;
  }

// ************** Initialize random seed (defined in library time.h) ***********
  srand(time(NULL));
  Float_t randvalue;

  ifstream inFileC;
  bool eof();

  bufana = fopen("ptest_bufana_31015_1016_1.dat","wt");   
       
// %%%%%%%%%%%%%%%%%%%  Here START the loop on RUN (FILES) %%%%%%%%%%%%%%%%%%%%%
 
for(loopf=nfile1;loopf<nfile2+1;loopf++){     

// %%%%%%%%%%%%%%%%%%%%%%%% Example input : many files %%%%%%%%%%%%%%%%%%%%%%%%%

    if(loopf==1){
      sprintf(infilnam,"DAT001000"); //write the name of data file which is generated by corsika
      num_ev_ene=0;  // Number of events of the same range of energy
      num_ev_ene_g30=0; //ak Num. ev. same energy range (gamma=-3.0 put knee)      
      loopf_ene=0;   //ak Range of energy
      enerefPeV=3; //ak Reference energy for this range in PeV (lowervalue)        
    }     
    if(loopf==2){
      sprintf(infilnam,"DAT001001");
      loopf_ene=0;   //ak Range of energy
      enerefPeV=3; //ak Reference energy for this range in PeV (lowervalue)        
    }  
    
     inFileC.open(infilnam,ios::binary); 
     fprintf (bufana,"\n %d) ====== OPEN FILE : %s ======== \n",loopf,infilnam);
     cout << loopf << ")  OPEN FILE : " << infilnam << endl;     
     nev = 0 ; // Event number given by Corsika buf[n][1]
     num_ev_run=0;  // Number of events and event number in the run (file)
     num_ev_test=0;
     num_ev_run_g30=0; //ak  Num. ev. with gamma=-3.0 in the run (put the knee)
          
  if (inFileC.bad()){
    fprintf (bufana,"Open input file ERROR unable to open it : loopf = %d \n",loopf);
    return(1);
  }

Int_t whilecount=0 ;

// %%% loop to read all the events in the open file DAT%%%%%% or input file %%% 

  while(1) {

    whilecount++;
    
// %%%%%%%%%%%%%%%%%%%% Here is the START to READ a record %%%%%%%%%%%%%%%%%%%%%

    for (int n=0;n<21;n++) { 
      for (int l=0;l<273;l++) {

        buf[n][l]=0;
   
	if(n==0 && l==0){

	  inFileC.read((char*)&firstlast,sizeof(Float_t));
	}

	inFileC.read((char*)&buf[n][l],sizeof(Float_t)); 

	if(n==20 && l==272){
	  inFileC.read((char*)&firstlast,sizeof(Float_t));
	}

        if(l==0)head.header = buf[n][l] ;
 
      }   
    }     

// *********************** Here is the END to READ a record *******************

// %%%%%%%%%%%%%%%%%%%%%%%%% START ANALYSING a record %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for (int n=0;n<21;n++) { 

      head.header = buf[n][0] ;
   
            
      if(strncmp(head.headname,"RUNH",4)==0) {        
        run = buf[n][1];
        slope = buf[n][15];
        emin = buf[n][16];
        emax = buf[n][17];
      }        // end RUNH


      if(strncmp(head.headname,"EVTH",4)==0) {
       
        if(nev==0) primev = buf[n][1];
	nev  = buf[n][1];
        num_ev_test ++;	

        eneprim  = buf[n][3];  // primary energy

 
        tagg30 = 0;  //ak
        eneprimPeV=eneprim/1000000;  //ak Primary Energy in PeV
        rapene=TMath::Power(eneprimPeV,deltakneegen)*TMath::Power(enerefPeV,deltagenknee); //ak ratio to choose
 
       // Take a random number between 0 and 1
        randvalue = rand()%32000+1; 
        randvalue /= 32000; 
	if(rapene>=randvalue){     //ak

          tagg30=1;                //ak
          num_ev_run_g30 ++;       //ak
          num_ev_ene_g30 ++; //ak Num.ev.in the same energy range g30

          if(num_ev_run_g30%100==0){  
            printf ("Loopf= %d, Loopf_ene= %d,  NUM_EV_RUN_G30= %d, nev = %f, NUM_SAME_ENE_G30= %d  Max.n.ev.same en.range= %d \n",loopf,loopf_ene,num_ev_run_g30,nev,num_ev_ene_g30,nevf_ene[loopf_ene]);
            printf ("nevf_ene[0]= %d, nevf_ene[1]= %d, nevf_ene[2]= %d, nevf_ene[3]= %d, nevf_ene[4]= %d \n",nevf_ene[0],nevf_ene[1],nevf_ene[2],nevf_ene[3],nevf_ene[4]);

           printf ("*******************************************\n \n");
	  }
        } // END  if(rapene>=randvalue)
//              End part after the knee   ak        
        
        
        if (nev==primev) {
 	  fprintf(bufana,"Event Number nev : %f \n",nev);
	  fprintf(bufana,"Run Number = %f, Slope en. spectrum = %f \n",run,slope);
	  fprintf(bufana,"Emin = %f, Emax = %f \n",emin,emax);
       	  fprintf(bufana,"Range of theta: %f %f\n",buf[n][80],buf[n][81]);
	  fprintf(bufana,"Range of phi  : %f %f\n",buf[n][82],buf[n][83]);
	  fprintf(bufana,"NKG radial distr. range in cm %f\n",buf[n][146]);  // 147 ?
          fprintf(bufana,"Event Number = %d, Prim. Particle = %f, Energy  %f \n",num_ev_run+1,part_id,eneprim);
        }

        // Take randomly the core coordinate xco,yco inside the sampling area
        randvalue = rand()%32000+1;  //integer random number between 1-32000 
        randvalue /= 32000; /*random value between ~0 and 1.*/
	xco = randvalue * ltot;
        randvalue = rand()%32000+1;  //integer random number between 1-32000 
        randvalue /= 32000; /*random value between ~0 and 1.*/
	yco = randvalue * ltot;


      	// Initialize at the beginning of the event the vectors for muons
        nmuali=0;
        nmuinali=0; 

	nmualix=0;
	nmuinalix=0;

	nmualiy=0;
	nmuinaliy=0;

        for(Int_t imu=0;imu<maxmu;imu++){
	  muplus[imu]=0;
	  muminus[imu]=0;
          tagmucell[imu]=0;
	  emucell[imu]=0;          
	  tetamucell[imu]=0;
          phimucell[imu]=0;
        }
                  

	// Initialize muon counters
        nmutot=nmut15=0;
        for (int k=0;k<LOOPCORE;k++) 
	  nmutdis[k] = 0;
	for (int k=0;k<ILOOP;k++) 
	  nmu[k] = tnmu[k] = emur[k] = temur[k] = 0.;

      }   // END  if(strncmp(head.headname,"EVTH",4)==0)  END EVTH

      if(strncmp(head.headname,"EVTE",4)==0) {
       
        sizenkg = buf[n][184]; 

        nmuev=int(nmuali);  
        nmuvec=nmuev-1;   
   	
	nmuevx=int(nmualix);
	nmuvecx=nmuevx-1;	

	nmuevy=int(nmualiy);
	nmuvecy=nmuevy-1;


//    CUT THE EVENT TO HAVE A POWER LAW GAMMA=gammaknee (tagg30=1) 
//             instead gamma=-2.7 that is the generation
      if(tagg30==1){   //ak
  	
        if(nmuev==0){ 
          nmuvec=0;     
          tagmu[nmuvec]=0;
          tetamuev[nmuvec]=-1000;
          phimuev[nmuvec]=-1000;
	  emuev[nmuvec]=-1000;
        } 

        if(nmuevx==0){ 
          nmuvecx=0;     
          tagmu[nmuvecx]=0;
	  tetamuev[nmuvecx]=-1000;
          phimuev[nmuvecx]=-1000;
        } 
       
        if(nmuevy==0){ 
          nmuvecy=0;     
          tagmu[nmuvecy]=0;
	  tetamuev[nmuvecy]=-1000;
          phimuev[nmuvecy]=-1000;
        }    
	
        if(nmuev==1){
          tagmu[nmuvec]=tagmucell[nmuvec];
          emuev[nmuvec]=emucell[nmuvec];
	  tetamuev[nmuvec]=tetamucell[nmuvec];
          phimuev[nmuvec]=phimucell[nmuvec];
        } 

        if(nmuevx==1){
          tagmu[nmuvecx]=tagmucell[nmuvecx];
          tetamuev[nmuvecx]=tetamucell[nmuvecx];
          phimuev[nmuvecx]=phimucell[nmuvecx];
        } 

        if(nmuevy==1){
          tagmu[nmuvecy]=tagmucell[nmuvecy];
          tetamuev[nmuvecy]=tetamucell[nmuvecy];
          phimuev[nmuvecy]=phimucell[nmuvecy];
        } 

	if(nmuev>1){ 
	 
          tagmu[0]=tagmucell[0];             
          emuev[0]=emucell[0]; 
	  tetamuev[0]=tetamucell[0];
          phimuev[0]=phimucell[0];
           
	  for(Int_t ias=1;ias<nmuev;ias++){
	   
            tagmu[ias]=tagmucell[ias];
            emuev[ias]=emucell[ias];
	    tetamuev[ias]=tetamucell[ias];
            phimuev[ias]=phimucell[ias];
       
	  }  
	} 

	if(nmuevx>1){ 
	 
          tagmu[0]=tagmucell[0];             
          tetamuev[0]=tetamucell[0];
          phimuev[0]=phimucell[0];
           
	  for(Int_t ias=1;ias<nmuevx;ias++){
	   
            tagmu[ias]=tagmucell[ias];
            tetamuev[ias]=tetamucell[ias];
            phimuev[ias]=phimucell[ias];
       
	  }  
	} 

	if(nmuevy>1){ 
	 
          tagmu[0]=tagmucell[0];             
          tetamuev[0]=tetamucell[0];
          phimuev[0]=phimucell[0];
           
	  for(Int_t ias=1;ias<nmuevy;ias++){
	   
            tagmu[ias]=tagmucell[ias];
            tetamuev[ias]=tetamucell[ias];
            phimuev[ias]=phimucell[ias];
       
	  }  
	} 

            if(nmuev>=0) ev_tree->Fill();
	    num_ev_run ++;

        if(num_ev_run%100==0||num_ev_run==nevf[loopf]){
          printf ("Loopf= %d, nev = %f, NUM_EV_RUN= %d  Total_events_run= %d \n",loopf,nev,num_ev_run,nevf[loopf]);
          printf ("*******************************************\n \n");
	}
               

	    // Check on tree variabiles 
            if(nmuev>0&&loopf==1&&num_ev_run<=100){
              fprintf(bufana,"********* Start Check on TREE variables ******\n ");
              fprintf(bufana,"Ev.number= %d, N.mu large area= %d   N.mu in Alice= %d \n",num_ev_run,nmuev,nmuinali);              

	      for(Int_t ias=0;ias<nmuev;ias++){ 
                fprintf(bufana,"Muon num.= %d, tagmu[ias]= %d \n",ias,tagmu[ias]);                                     
              } 
            } // END if(nmuev>0)
        } //ak END if(tagg30==1)  do not write the event if tagg30=0    
      }   // END if(strncmp(head.headname,"EVTE",4)==0)   END EVTE event end


      if(strncmp(head.headname,"RUNE",4)==0) {
	
	fprintf(bufana," ============= Found RUNE ============ Block \n");

      	fprintf(bufana,"N. events in the run = %d, Total N. events before this run = %d\n",num_ev_run,num_ev_tot);

	fprintf(bufana,"==============END RUNE=================== \n");
       
      }   //  END  if(strncmp(head.headname,"RUNE",4)==0)   end RUNE  run end
    
// %%%%%%%%%%%%%%%%% Analysis of the particles function AnaEv %%%%%%%%%%%%%%%%%%

      if(strncmp(head.headname,"RUNH",4)!=0 && strncmp(head.headname,"EVTH",4)!=0 &&
	 strncmp(head.headname,"EVTE",4)!=0 && strncmp(head.headname,"RUNE",4)!=0)
	{AnaEv(n); /* n is the block */
	} 
    
    } // END for (int n=0;n<21;n++) End Analazing a record buf[21][273]   


    if(num_ev_ene_g30>=nevf_ene[loopf_ene]+1||inFileC.eof()){  //ak

      num_ev_tot += num_ev_run;   //ak total number of events
//      num_ev_ene += num_ev_run;   // total num. of ev. in the same energy range
      num_ev_tot_g30 += num_ev_run_g30;   //ak total number of events g30
//      num_ev_ene_g30 += num_ev_run_g30;   // tot.num.of ev.same en.range g30

      fprintf(bufana," ================ SUMMARY RUN ======================== \n");
      fprintf(bufana,"Loopf  = %d , Num.ev.run g30 = %d , Num.ev.run = %d \n",loopf,num_ev_run_g30,num_ev_run);
      fprintf(bufana,"Loopf_ene = %d, MAX ev.ene range g30 = %d , NUM.EV.SAME ENERGY RANGE g30 = %d, Num.ev.consumed same ene = %d \n",loopf_ene,nevf_ene[loopf_ene],num_ev_ene_g30,num_ev_ene);
     fprintf(bufana," Tot.ev consumed = %d, TOT.EV g30 = %d \n",num_ev_tot,num_ev_tot_g30); 
      fprintf(bufana," ================ END SUMMARY RUN ======================== \n");
     
      cout << "Loop = " <<loopf << "   Num. ev run g30 = " << num_ev_run_g30 << "  Num. ev. run = " << num_ev_run << "   Tot. ev consumed = " << num_ev_tot << "  TOTAL EVENTS gammaknee = " << num_ev_tot_g30 << endl;  
                               
      break;     // stop loop on this run (file)
    
    }     // END if(num_ev_run==nevf[loopf]||inFileC.eof())

  }    // end while(1) loop to read all the events in the open file 

  printf ("End Read Events\n");
  inFileC.close();

  }   // end loop on input files : loopf

   ftree1->cd();
  ftree1->Write();
  ftree1->Close();
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%% ANALYSIS OF MUONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void AnaEv(int nb)
{
  int j,jj;
  int ip,idp;
  Float_t part[word_part][maxpart_db];
  Float_t pxmu,pymu,pzmu,xmu,ymu,zmu;
  Float_t enemu,mommu,costetamu,tetamu,phimu;
  Float_t ellemu,emmemu,ennemu,tvar;
  Float_t xmualilevel,ymualilevel;
  Int_t nmuv,nmuvx,nmuvy;
  Int_t chmuon,chmuonp,chmuonm;

 //   printf("Enter in AnaEv \n");

 // maxpart = 39     word_part = 7   39x7=273 sub-block

  for (j=0;j<maxpart_db;j++) {
    for (jj=0;jj<word_part;jj++) {

      ip = j*7 + jj ; // j is the particle, jj is the word 
      // printf("n=%d l=%d buf=%f \n",nb,ip,buf[nb][ip]);  
      part[jj][j] = buf[nb][ip];   // particle sub-block
    }                                   
   
    idp = (int)(part[0][j]/1000); // particle identifier
     
 // printf("Particle id = %d \n",idp);
 // if (idp!=5 && idp!=6 && idp!=75 && idp!=76)  //non muons

    if (idp!=0) {               // in case of last block not completly filled   
      if(idp==5 || idp==6) {                   
                              // printf("Inside muons Particle id = %d \n",idp);
	pxmu = part[1][j];
	pymu = part[2][j];
	pzmu = part[3][j];
	xmu  = part[4][j];//Corsika x coord.mu at surf.(xcore=0,ycore=0,zcore=0) 
	ymu  = part[5][j];//Corsika y coord.mu at surf.(xcore=0,ycore=0,zcore=0)
	zmu  = 0.;        //Put surface at z=0    
	enemu = sqrt(pxmu*pxmu + pymu*pymu + pzmu*pzmu + MASSMU*MASSMU);
	mommu = sqrt(pxmu*pxmu + pymu*pymu + pzmu*pzmu);  
	costetamu = pzmu/mommu;
	tetamu = acos(costetamu);
	phimu = atan(pymu/pxmu);

        if(idp==5) chmuon=1;
	if(idp==6) chmuon=-1;

	if (pxmu<0.&&pymu>=0.) 
	  phimu=phimu+PIGR;     
	if (pxmu<0.&&pymu<0.)
	  phimu = phimu+PIGR;   
	if (pxmu>=0. && pymu<0.)
	  phimu = phimu+2*PIGR; 

	if (enemu>0.) {

	  nmutot++;                           

	  if (enemu>1) {       
	  
	    nmut15++;

//  Muon coordinates at Alice level

            ellemu = sin(tetamu)*cos(phimu);
	    emmemu = sin(tetamu)*sin(phimu);
	    ennemu = cos(tetamu);
	    tvar = (zdeep-zmu)/ennemu;

// xmualilevel ymualilevel coord. of muon at Alice level (40 m underground) 
// keeping into account the core of the shower xco,yco
 
	    xmualilevel = xmu + ellemu * tvar + xco;
	    ymualilevel = ymu + emmemu * tvar + yco;

// Check if muon is inside an area 6 x 6 m**2 centered in Alice and write the variables
              if(xmualilevel>=xylowlargearea&&xmualilevel<=xyhighlargearea){
              if(ymualilevel>=xylowlargearea&&ymualilevel<=xyhighlargearea){
              nmuali+=1;    
              nmuv=int(nmuali)-1;	      
              emucell[nmuv]=enemu; 
	      tetamucell[nmuv]=tetamu*57.297;
	      phimucell[nmuv]=phimu*57.297;
	      idchmu[nmuv]=chmuon;
//   Tag mu. If muons inside an area 5 x 3 m**2 and E>16/(cos(teta)
//                 tagmu[nmuv]=1 otherwise tagmu[nmuv]=0	    
              tagmucell[nmuv]=0; 
              if(ymualilevel>=(yalilow3m)&&ymualilevel<=(yalihigh3m)){
                if(xmualilevel>=(xalilow)&&xmualilevel<=(xalihigh)){
                  if (enemu>(emucut/costetamu)) {
                    tagmucell[nmuv]=1;
                    nmuinali+=1; 
             	  }	
                }
               }
             } // END if(ymualilevel>=(yalilow3m-150).....
            } // END if(xmualilevel>=(xalilow-50) .....  
	  }   
	}   // end if (enemu>0.)   
      }    // end if mu+ or mu-    
 
     if(idp==5) {                   

	chmuonp=1;
	
       	if (enemu>0.) {

	  nmutot++;                           

	  if (enemu>1) {       
	  
	      nmut15++;
	      if(xmualilevel>=xylowlargearea&&xmualilevel<=xyhighlargearea){
              if(ymualilevel>=xylowlargearea&&ymualilevel<=xyhighlargearea){
              nmualix+=1;    
              nmuvx=int(nmualix)-1;	      
	      muplus[nmuvx]=chmuonp;   
	       }  
	     }   
	  }   
	}  //end if (enemu>0.) 
     }   //end if (idp==5)
    
    if(idp==6) {                   

       chmuonm=-1;
	
       	if (enemu>0.) {

	  nmutot++;                           

	  if (enemu>1) {       
	     
              nmut15++;
              if(xmualilevel>=xylowlargearea&&xmualilevel<=xyhighlargearea){
              if(ymualilevel>=xylowlargearea&&ymualilevel<=xyhighlargearea){
              nmualiy+=1;    
              nmuvy=int(nmualiy)-1;	      
	      muminus[nmuvy]=chmuonm;   
	     } // END if(ymualilevel>=(yalilow3m-150)..... 
	    } // END if(xmualilevel>=(xalilow-50) .....
 
	 }  //end if mu above threshold   
       }   //end if (enemu>0.) 
     }    //end if (idp==6)
    }    //end if idp!=0
  }     //end for j=0,maxpart_db
}      // AnaEv(int nb)       




