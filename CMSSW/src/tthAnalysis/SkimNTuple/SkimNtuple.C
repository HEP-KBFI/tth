#define SkimNtuple_cxx
#include "SkimNtuple.h"
#include <iostream>
#include <cmath>
#include <algorithm>    // std::find
#include <vector>       // std::vector
#include <TH1.h>   
#include <TLorentzVector.h>

using namespace std;

void SkimNtuple::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L SkimNtuple.C
//      Root > SkimNtuple t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


   if (fChain == 0) return;

   nentries = fChain->GetEntries();
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) { //evt
   //for (Long64_t jentry=0; jentry<20000;jentry++) { //evt
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


     ////////////////////////////////////////
     //Make a skim w/ at leat 1 tight lepton
     ///////////////////////////////////////

      nbmu=0; nbe=0; nbtau=0;

      for ( int i=0; i<np; i++){
	      if(abs(pid[i])==13 && (id[i]&idtmu) ==idtmu  ) nbmu++;
	      if(abs(pid[i])==11 && (id[i]&idte)  ==idte   ) nbe++;
	      if(abs(pid[i])==15 && (id[i]&idltau)==idltau ) nbtau++;
      }
      nev++;

      if( nbmu>0) ntmu++; if( nbe>0) nte++; if( nbtau>1) ntau++;

      if( (nbmu>0 || nbe>0) && nbtau>0 ) skimtree->Fill();
      //if( nbmu>0 ) skimtree->Fill();


     ///////////////////////////////////
     //GenLevel study for tau eff.
     ///////////////////////////////////

/*
      //With this method I got 50% eff for a loose iso cut
      //according to the expert Tau::genJet() match ONLY hadronic tau to the generator
      for(int i=0; i<np; i++){

	      if( ! (abs(pid[i])==15 && (id[i]&1<<10)==1<<10)) continue;

	      rtau.SetPxPyPzE(px[i], py[i], pz[i], e[i]);

	      if(!(rtau.Pt()>20 && fabs(rtau.Eta())< 2.3)) continue;

	      h0->Fill(rtau.Pt());

	      if((id[i]&1<<0)==1<<0) h1->Fill(rtau.Pt());

	      if((id[i]&1<<1)==1<<1) h2->Fill(rtau.Pt());

	      if((id[i]&1<<2)==1<<2) h3->Fill(rtau.Pt());

	      if((id[i]&1<<3)==1<<3) h4->Fill(rtau.Pt());

	      if((id[i]&1<<4)==1<<4) h5->Fill(rtau.Pt());

	      if((id[i]&1<<8)==1<<8) h6->Fill(rtau.Pt());
      }
*/

/*
      //with this method I also got 50% eff for a loose cut
      //get all hadronic gen taus index
      vgtau.clear();  
      for ( int i=0; i<gnp; i++){ if ( abs(gpid[i])==15 && gstatus[i]==2 ) vgtau.push_back(i); }

      //remove lepton taus from all taus and keep hadronic taus
      for ( int i=0; i<gnp; i++){
	      if( abs(gpid[i])==11 || abs(gpid[i])==13 ){
		      igtau = find(vgtau.begin(), vgtau.end(), gmother[i]) ; // check if the lepton mother is a tau (include in vgtau)
		      if ( igtau!=vgtau.end() )  vgtau.erase(igtau); // if yes, remove this tau in vgtau
	      }
      }

      //Get the reco hadronic tau index
      vrtau.clear(); 
      for (int i=0; i<np; i++){ if( abs(pid[i])==15) vrtau.push_back(i);}

      //recoTau+GenTau matching

	
      for(int j : vgtau){//gen

	      Lgtau.SetPxPyPzE( gpx[j], gpy[j], gpz[j], ge[j] );
	      if( !(Lgtau.Pt()>20 && fabs(Lgtau.Eta())<2.3 )) continue;

	      for(int i : vrtau){//reco

		      Lrtau.SetPxPyPzE( px[i], py[i], pz[i], e[i] );
		      if( !(Lrtau.Pt()>20 && fabs(Lrtau.Eta())<2.3 )) continue;

		      if( !(Lrtau.DeltaR(Lgtau)<0.2 && (Lrtau.Pt()-Lgtau.Pt())<5)) continue; //recoMathGen
	    
		      h0->Fill(Lrtau.Pt());
		      if( (id[i]&1<<0)==1<<0) h1->Fill(Lrtau.Pt());
		      if( (id[i]&1<<1)==1<<1) h2->Fill(Lrtau.Pt());
		      if( (id[i]&1<<2)==1<<2) h3->Fill(Lrtau.Pt());
		      if( (id[i]&1<<3)==1<<3) h4->Fill(Lrtau.Pt());
	      }

      }
*/


   }//evt

   cout << "nentries: " <<nentries<< ", nev,: "<<nev <<endl;

}

