#define AnalyseNtuple_cxx
#include "AnalyseNtuple.h"
#include <iostream>
#include <cmath>

using namespace std;

void AnalyseNtuple::Loop()
{
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry

    nentries = fChain->GetEntries();
    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) { //evt
    //for (Long64_t jentry=0; jentry<20;jentry++) { //evt

        // check that the trees are setup properly
	Long64_t ientry = LoadTree(jentry); if (ientry < 0) break; 
	nb = fChain->GetEntry(jentry);   nbytes += nb; if (Cut(ientry) < 0) continue;

        // clean particle's vector while a new event begin
	std::vector< vector<int> > cleanv {vtau, vtmu, vlmu, vvmu, vte, vle, vve, vj}; for( auto i : cleanv) i.clear(); 
        // fill particle's vector depending on their ID (ele=11, muon=13, tau=15, jet=81) and passing cuts (loose, medium , tight) 
	for ( int i=0; i <np; i++){
	    if( fabs(pid[i])==11 ){ 
		if( (id[i]&idte)==idte && disc[i]>0.5 ) vte.push_back(i);	      
		if( (id[i]&idle)==idle && (id[i]&1<<0)!=1<<0 && disc[i]>0.5) vle.push_back(i);
		if( (id[i]&idve)==idve && (id[i]&1<<0)!=1<<0 && (id[i]&1<<1)!=1<<1 ) vve.push_back(i); //disc to be study!!
	    }
	    if( abs(pid[i])==13 ){
		if( (id[i]&idtmu)==idtmu )  vtmu.push_back(i); 
		if( (id[i]&idlmu)==idlmu && (id[i]&1<<0)!=1<<0 ) vlmu.push_back(i);
		if( (id[i]&idvmu)==idvmu && (id[i]&1<<0)!=1<<0 && (id[i]&1<<1)!=1<<1  ) vvmu.push_back(i); 
	    }
	    if( abs(pid[i])==15 && (id[i]&dlll)==dlll ) vtau.push_back(i);
	    if(pid[i]==81 && (id[i]&idlj)==idlj ) vj.push_back(i);
	}

	nev++; 

	// pileup rewighting for MC only
	if( sID==1<<14 || sID==1<<15 ) wpu=1.; else wpu=LumiWeights_.weight( float(npu) );
	//total weight to apply for each event
	w=wpu;        


	//////////////////////////
	// single muon analysis //
	//////////////////////////

	if( trig==1) { 
	    if( vtmu.size()<1 ) continue;
	    Ltmu.SetPxPyPzE(px[vtmu[0]], py[vtmu[0]], pz[vtmu[0]], e[vtmu[0]]);
	    if( Ltmu.Pt()<26 ) continue;
	    if( iso[vtmu[0]] >0.12 ) continue;
	    if( vtau.size() <2 ) continue;

	    metpt=sqrt(metx*metx+mety*mety); 
	    //if( metpt < 40 ) continue;

	    //qcd selection
	    ra=true; rb=false; rc=false; rd=false;
	    //ra = metpt>40. && iso[vtmu[0]]<0.1 ;    //signal region A
	    rb = metpt<40. && iso[vtmu[0]]<0.12 ;      //qcd region B
	    rc = metpt<40. && iso[vtmu[0]]>0.12 ;      //qcd region C
	    rd = metpt>40. && iso[vtmu[0]]>0.12 ;    //signal region D
	    region.clear(); region.push_back(ra); region.push_back(rb); region.push_back(rc); region.push_back(rd); 

	    // tau removal form jet collection
	    vjclean.clear(); vbclean.clear();
	    for (int i : vj ){  

		Lj.SetPxPyPzE(px[i], py[i], pz[i], e[i]);

		findtau=false;
		for (int j : vtau ){ Ltau.SetPxPyPzE(px[j], py[j], pz[j], e[j]); if( Lj.DeltaR(Ltau)<0.3 ) findtau=true; } 

		if ( findtau) continue; 

		if ( fabs(Lj.Eta()) >2.5 ) continue; 

		for (int k=0; k<region.size(); k++) if(region[k]) {
		    hjpt_sm[k]->Fill(Lj.Pt(), w); hjeta_sm[k]->Fill(Lj.Eta(), w); hzeta_sm[k]->Fill(disc[i], w);
		}

		if     ( Lj.Pt()>40 && disc[i]<0.679 ) vjclean.push_back(i);

		else if( Lj.Pt()>30 && disc[i]>0.679 ) vbclean.push_back(i); 

	    }
	    //if ( vjclean.size() <1 ) continue;
	    //if ( vbclean.size() <1 ) continue;

	    //Lj0.SetPxPyPzE(px[vjclean[0]], py[vjclean[0]], pz[vjclean[0]], e[vjclean[0]]);
	    //Lj1.SetPxPyPzE(px[vjclean[1]], py[vjclean[1]], pz[vjclean[1]], e[vjclean[1]]);
	    //Lb0.SetPxPyPzE(px[vbclean[0]], py[vbclean[0]], pz[vbclean[0]], e[vbclean[0]]);
	    //Lb1.SetPxPyPzE(px[vbclean[1]], py[vbclean[1]], pz[vbclean[1]], e[vbclean[1]]);
	    Ltau0.SetPxPyPzE(px[vtau[0]], py[vtau[0]], pz[vtau[0]], e[vtau[0]]);

	    for (int i=0; i<region.size(); i++) if(region[i]) {
		hnvtx_sm[i]   ->Fill( nvtx, w);
		hmetpt_sm[i]  ->Fill( metpt, w );
		hmupt_sm[i]   ->Fill( Ltmu.Pt(), w );     
		hmueta_sm[i]  ->Fill( Ltmu.Eta(), w ); 
		hmuphi_sm[i]  ->Fill( Ltmu.Phi(), w ); 
		hjmult_sm[i]  ->Fill( vjclean.size(), w ); 
		hbmult_sm[i]  ->Fill( vbclean.size(), w ); 
		htau0pt_sm[i] ->Fill( Ltau0.Pt(), w );  
		htau0eta_sm[i]->Fill( Ltau0.Eta(), w );
		hdrtaumu_sm[i] ->Fill(Ltau0.DeltaR(Ltmu), w );     
	    }


	}
	//hj0pt_sm[0]  ->Fill(Lj0.Pt(), w ); 
	//hj0eta_sm[0] ->Fill(Lj0.Eta(), w ); 
	//hj0e_sm[0]   ->Fill(Lj0.E(), w ); 
	//hj1pt_sm[0]  ->Fill(Lj1.Pt(), w ); 
	//hj1eta_sm[0] ->Fill(Lj1.Eta(), w ); 
	//hj1e_sm[0]   ->Fill(Lj1.E(), w );
	//hjdr01_sm[0] ->Fill(Lj0.DeltaR(Lj1), w );     
	//h2jm_sm[0]   ->Fill((Lj0+Lj1).M(), w );
	/*
	   if( vtau.size()>1){
	   Ltau1.SetPxPyPzE(px[vtau[1]], py[vtau[1]], pz[vtau[1]], e[vtau[1]]);
	   htau1pt_sm[0] ->Fill( Ltau1.Pt(), w );   
	   htau1eta_sm[0]->Fill( Ltau1.Eta(), w ); 
	   hdr2tau_sm[0] ->Fill(Ltau0.DeltaR(Ltau1), w );
	   }
	 */
	// }

	/*
	//el
	//if( trig==2 && vte.size()==1 && vle.size()==0 && vtmu.size()==1 && vlmu.size()==0){
	if( trig==2 && vte.size()>0 ){
	Lte.SetPxPyPzE(px[vte[0]], py[vte[0]], pz[vte[0]], e[vte[0]]);
	if( Lte.Pt()<30 ) continue;
	if( iso[vte[0]] <0.10 ) continue;

	hmetpt_se[0]   ->Fill(metpt, w );
	hnvtx_se[0]    ->Fill(nvtx, w);
	hept_se[0]     ->Fill(Ltmu.Pt(), w );     
	heeta_se[0 ]   ->Fill(Ltmu.Eta(), w ); 
	hjmult_se[0]   ->Fill(vjclean.size(), w ); 
	hbmult_se[0]   ->Fill(vbclean.size(), w ); 
	hj0pt_se[0]    ->Fill(Lj0.Pt(), w ); 
	hj0eta_se[0]   ->Fill(Lj0.Eta(), w ); 
	hj0e_se[0]     ->Fill(Lj0.E(), w ); 
	hj1pt_se[0]    ->Fill(Lj1.Pt(), w ); 
	hj1eta_se[0]   ->Fill(Lj1.Eta(), w ); 
	hj1e_se[0]     ->Fill(Lj1.E(), w );
	hjdr01_se[0]   ->Fill(Lj0.DeltaR(Lj1), w );     
	h2jm_se[0]     ->Fill((Lj0+Lj1).M(), w );
	htau0pt_se[0]  ->Fill( Ltau0.Pt(), w );  
	htau0eta_se[0] ->Fill( Ltau0.Eta(), w );

	if( vtau.size()>1){
	Ltau1.SetPxPyPzE(px[vtau[1]], py[vtau[1]], pz[vtau[1]], e[vtau[1]]);
	htau1pt_se[0] ->Fill( Ltau1.Pt(), w );   
	htau1eta_se[0]->Fill( Ltau1.Eta(), w ); 
	hdr2tau_se[0] ->Fill(Ltau0.DeltaR(Ltau1), w );
	}
	}
	 */



	/*

	///////////////////////////////////
	// l+j (1muon, 2had taus) Analysis   
	///////////////////////////////////   

	if ( trig==1 ){ nsmutrig++;              
	if( vtmu.size()!=1 ) continue ; n1tmu++;
	if( vlmu.size()!=0 ) continue;  nolmu++;
	if( vle.size()!=0)   continue;  noel++;

	if(!(Ltmu.Pt()>30 && fabs(Ltmu.Eta())<2.1 && iso[vtmu[0]]<0.1 )) continue; ntmu++;
	if( vtau.size()<2)     continue; n2tau++;

	hSLmupt->Fill(Ltmu.Pt());     
	hSLmueta->Fill(Ltmu.Eta()); 
	htau0pt->Fill( Ltau0.Pt());  
	htau0eta->Fill( Ltau0.Eta());
	if( vtau.size()>2){
	Ltau1.SetPxPyPzE(px[vtau[1]], py[vtau[1]], pz[vtau[1]], e[vtau[1]]);
	htau1pt->Fill( Ltau1.Pt());   htau1eta->Fill( Ltau1.Eta()); 
	hallpt ->Fill((Ltmu+Lj0+Lj1+Lb0+Lb1+Ltau0+Ltau1+Lmet).Pt()); //pt(muon, jet, tau, met)
	hallm  ->Fill((Ltmu+Lj0+Lj1+Lb0+Lb1+Ltau0+Ltau1+Lmet).M()); //mass(muon, jet, tau, met)
	}


	hjmult ->Fill(vjclean.size()); 
	h2jm   ->Fill((Lj0+Lj1).M());

	hb0pt  ->Fill(Lb0.Pt()); hb0e->Fill(Lb0.E()); //hbzeta->Fill(disc[i]);
	hb1pt  ->Fill(Lb1.Pt()); hb1e->Fill(Lb1.E()); //hbzeta->Fill(disc[i]);
	hbmult ->Fill(vbclean.size());
	h2bm   ->Fill((Lb0+Lb1).M());
	}


	/////////////////////////////////////////////
	// l+j (1muon, 2taus (1had + 1muon) Analysis   
	/////////////////////////////////////////////   

	if ( trig==1 ){ 
	if( vtmu.size()!=2 ) continue; 
	if( vlmu.size()!=0 ) continue; 
	if( vle.size()!=0)   continue; 

	m2_Ltmu1.SetPxPyPzE(px[vtmu[1]], py[vtmu[1]], pz[vtmu[1]], e[vtmu[1]]);
	if(!(Ltmu.Pt()>30 && fabs(Ltmu.Eta())<2.1 && iso[vtmu[0]]<0.1 )) continue; 
	if(!(m2_Ltmu1.Pt()>30 && fabs(m2_Ltmu1.Eta())<2.1 && iso[vtmu[1]]<0.1 )) continue;
	if( vtau.size()<1)     continue; 

	m2_hSLmu0pt  ->Fill( Ltmu.Pt());     
	m2_hSLmu0eta ->Fill( Ltmu.Eta()); 
	m2_hSLmu1pt  ->Fill( m2_Ltmu1.Pt());     
	m2_hSLmu1eta ->Fill( m2_Ltmu1.Eta()); 
	m2_htau0pt   ->Fill( Ltau0.Pt());   
	m2_htau0eta  ->Fill( Ltau0.Eta());

	m2_hmetpt ->Fill(metpt);

	m2_hj0pt  ->Fill(Lj0.Pt());   
	m2_hj0e   ->Fill(Lj0.E()); //hjzeta->Fill(disc[i]);
	m2_hj1pt  ->Fill(Lj1.Pt());   
	m2_hj1e   ->Fill(Lj1.E()); //hjzeta->Fill(disc[i]);
	m2_hjmult ->Fill(vjclean.size()); 
	m2_h2jm   ->Fill((Lj0+Lj1).M());

	m2_hb0pt  ->Fill(Lb0.Pt());   
	m2_hb0e   ->Fill(Lb0.E()); //hbzeta->Fill(disc[i]);
	m2_hb1pt  ->Fill(Lb1.Pt());  
	m2_hb1e   ->Fill(Lb1.E()); //hbzeta->Fill(disc[i]);
	m2_hbmult ->Fill(vbclean.size());
	m2_h2bm   ->Fill((Lb0+Lb1).M());
    }
    */

	/*
	//if is singleEl!!
	//njet>=2, tighte=1, veto
	if( vjclean.size()>=2 && vtau.size()>1 && vte.size()==1 && vle.size()==0 && vtmu.size()==0 && vlmu.size()==0){ //vvpt shoud be >20 GeV, here >10 GeV !!
	Lte.SetPxPyPzE(px[vte[0]], py[vte[0]], pz[vte[0]], e[vte[0]]);

	if(!(Lte.Pt()>30 && fabs(Lte.Eta())<2.5 && (fabs(Lte.Eta())<1.442 || fabs(Lte.Eta())>1.566) && iso[vte[0]]<0.1)) continue;

	hSLept->Fill(Lte.Pt());

	nte++;
	//isej=true;
	}
	 */

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////
	//    ll Analysis   //
	///////////////////////   
	//is2tmu=false; is2te=false; istmue=false; istlmu=false; istle=false; istmule=false; istelmu=false;
	/* 
	   pass0t=false; pass0l=false; pass1t=false; pass1l=false; 

	//di-tmu
	if( vtmue.size()>1 && vtee.size()==0 ){ 

	vl0.SetPxPyPzE(px[vtmu[0]], py[vtmu[0]], pz[vtmu[0]], e[vtmu[0]]);
	vl1.SetPxPyPzE(px[vtmu[1]], py[vtmu[1]], pz[vtmu[1]], e[vtmu[1]]);

	pass0t =vl0.Pt()>20 && fabs(vl0.Eta())<2.1 && vtmuiso[0] <0.1;
	pass0l =vl0.Pt()>10 && fabs(vl0.Eta())<2.3 && vtmuiso[0] <0.2;
	pass1t =vl1.Pt()>20 && fabs(vl1.Eta())<2.1 && vtmuiso[0] <0.1;
	pass1l =vl1.Pt()>10 && fabs(vl1.Eta())<2.3 && vtmuiso[0] <0.2;

	if( (pass0t&&pass1l || pass0l&&pass1t ) && pid[vtmu[0]]*pid[vtmu[1]]<0 ){

	//cout << "q1: "<< vtmuq[0] <<", q2: "<< vtmuq[1] << ", q1q2<0? "<< vtmuq[0]*vtmuq[1] << endl;
	n2tmu++;
	is2tmu=true;
	}
	}
	 */

	/*

	//di-tele
	//else if(( vtee.size()==2 ) && ( vlee.size()==0 && vtmue.size()==0 && vlmue.size()==0) ){
	else if( vtee.size()==2 && vtmue.size()==0 ){

	vl0.SetPxPyPzE(vtepx[0], vtepy[0], vtepz[0], vtee[0]);
	vl1.SetPxPyPzE(vtepx[1], vtepy[1], vtepz[1], vtee[1]);

	pass0t =vl0.Pt()>20 && fabs(vl0.Eta())<2.5 && (fabs(vl0.Eta())<1.442 || fabs(vl0.Eta())>1.566) && vteiso[0]<0.1; 
	pass0l =vl0.Pt()>10 && fabs(vl0.Eta())<2.5 && (fabs(vl0.Eta())<1.442 || fabs(vl0.Eta())>1.566) && vteiso[0]<0.2;
	pass1t =vl1.Pt()>20 && fabs(vl1.Eta())<2.5 && (fabs(vl1.Eta())<1.442 || fabs(vl1.Eta())>1.566) && vteiso[1]<0.1;
	pass1l =vl1.Pt()>10 && fabs(vl1.Eta())<2.5 && (fabs(vl1.Eta())<1.442 || fabs(vl1.Eta())>1.566) && vteiso[1]<0.2;

	if( ( pass0t&&pass1l || pass0l&&pass1t ) && vteq[0]*vteq[1]<0){

	n2te++;
	is2te=true;
	}
	}

	//tmu-tele
	//else if( (vtee.size()==1 && vtmue.size()==1) && (vlee.size()==0 && vlmue.size()==0) ){
	else if( vtee.size()==1 && vtmue.size()==1 ){

	vl0.SetPxPyPzE(vtepx[0], vtepy[0], vtepz[0], vtee[0]);
	vl1.SetPxPyPzE(vtmupx[0], vtmupy[0], vtmupz[0], vtmue[0]);

	pass0t =vl0.Pt()>20 && fabs(vl0.Eta())<2.5 && (fabs(vl0.Eta())<1.442 || fabs(vl0.Eta())>1.566) && vteiso[0]<0.1;
	pass0l =vl0.Pt()>10 && fabs(vl0.Eta())<2.5 && (fabs(vl0.Eta())<1.442 || fabs(vl0.Eta())>1.566) && vteiso[0]<0.2;
	pass1t =vl1.Pt()>20 && fabs(vl1.Eta())<2.1 && vtmuiso[0] <0.1;
	pass1l =vl1.Pt()>10 && fabs(vl1.Eta())<2.3 && vtmuiso[0] <0.2;

	if( (pass0t&&pass1l ||pass0l&&pass1t) && vteq[0]*vtmuq[0]<0 ){

	n2lep++;
	istmue=true;
	}
	}

	//t+l mu
	//else if( ( vtmue.size()==1 && vlmue.size()==1 ) && (vtee.size()==0 && vlee.size()==0) ){
	else if( vtmue.size()==1 && vlmue.size()==1 ){

	vl0.SetPxPyPzE(vtmupx[0], vtmupy[0], vtmupz[0], vtmue[0]);
	vl1.SetPxPyPzE(vlmupx[0], vlmupy[0], vlmupz[0], vlmue[0]);

	pass0t =vl0.Pt()>20 && fabs(vl0.Eta())<2.1 && vtmuiso[0] <0.1;
	pass0l =vl0.Pt()>10 && fabs(vl0.Eta())<2.3 && vtmuiso[0] <0.2;
	pass1t =vl1.Pt()>20 && fabs(vl1.Eta())<2.1 && vlmuiso[0] <0.1;
	pass1l =vl1.Pt()>10 && fabs(vl1.Eta())<2.3 && vlmuiso[0] <0.2;

	if( (pass0t&&pass1l || pass0l&&pass1t ) && vtmuq[0]*vlmuq[0]<0 ){

	n2mu++;
	istlmu=true;
	}
	}

	//t+l e
	//else if( ( vtee.size()==1 && vlee.size()==1 )  && (vtmue.size()==0 && vlmue.size()==0) ){
	else if( vtee.size()==1 && vlee.size()==1 ){

	vl0.SetPxPyPzE(vtepx[0], vtepy[0], vtepz[0], vtee[0]);
	vl1.SetPxPyPzE(vlepx[0], vlepy[0], vlepz[0], vlee[0]);

	pass0t =vl0.Pt()>20 && fabs(vl0.Eta())<2.5 && (fabs(vl0.Eta())<1.442 || fabs(vl0.Eta())>1.566) && vteiso[0]<0.1 ;
	pass0l =vl0.Pt()>10 && fabs(vl0.Eta())<2.5 && (fabs(vl0.Eta())<1.442 || fabs(vl0.Eta())>1.566) && vteiso[0]<0.2;
	pass1t =vl1.Pt()>20 && fabs(vl1.Eta())<2.5 && (fabs(vl1.Eta())<1.442 || fabs(vl1.Eta())>1.566) && vleiso[1]<0.1;
	pass1l =vl1.Pt()>10 && fabs(vl1.Eta())<2.5 && (fabs(vl1.Eta())<1.442 || fabs(vl1.Eta())>1.566) && vleiso[1]<0.2;

	if( ( pass0t&&pass1l || pass0l&&pass1t ) && vteq[0]*vleq[0]<0 ){

	    n2e++;
	    istle=true;
	}
}

//tight(mu) loose(e)
else if( vtmue.size()==1 && vlee.size()==1   && ( vlmue.size()==0 && vtee.size()==0) ){

    vl0.SetPxPyPzE(vtmupx[0], vtmupy[0], vtmupz[0], vtmue[0]);
    vl1.SetPxPyPzE(vlepx[0], vlepy[0], vlepz[0], vlee[0]);

    pass0t =vl0.Pt()>20 && fabs(vl0.Eta())< 2.1 && vtmuiso[0] <0.1;
    pass0l =vl0.Pt()>10 && fabs(vl0.Eta())< 2.3 && vtmuiso[0] <0.2;
    pass1t =vl1.Pt()>20 && fabs(vl1.Eta())<2.5 && (fabs(vl1.Eta())<1.442 || fabs(vl1.Eta())>1.566) && vteiso[1]<0.1;
    pass1l =vl1.Pt()>10 && fabs(vl1.Eta())<2.5 && (fabs(vl1.Eta())<1.442 || fabs(vl1.Eta())>1.566) && vteiso[1]<0.2;

    if( (pass0t&&pass1l || pass0l&&pass1t ) && vtmuq[0]*vleq[0]<0 ){

	nltmue1++;
	istmule=true;
    }
}

//tight(e) loose(mu)
else if( (vtee.size()==1  && vlmue.size()==1) && ( vlee.size()==0  && vtmue.size()==0) ){

    vl0.SetPxPyPzE(vtepx[0], vtepy[0], vtepz[0], vtee[0]);
    vl1.SetPxPyPzE(vlmupx[0], vlmupy[0], vlmupz[0], vlmue[0]);

    pass0t =vl0.Pt()>20 && fabs(vl0.Eta())<2.5 && (fabs(vl0.Eta())<1.442 || fabs(vl0.Eta())>1.566) && vteiso[0]<0.1;
    pass1l =vl1.Pt()>10 && fabs(vl1.Eta())<2.3 && vlmuiso[0] <0.2;  
    pass0l =vl0.Pt()>10 && fabs(vl0.Eta())<2.5 && (fabs(vl0.Eta())<1.442 || fabs(vl0.Eta())>1.566) && vteiso[0]<0.2;
    pass1t =vl1.Pt()>20 && fabs(vl1.Eta())<2.1 && vlmuiso[0] <0.1; 

    if(  ( pass0t&&pass1l || pass0l&&pass1t ) && vteq[0]*vlmuq[0]<0 ){ 

	nltmue2++;
	istelmu=true;
    }
}

//other
else ntrash++;
*/

/////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////
// event categories: SL //
//////////////////////////
/*
//1t e/mu 

if( (is2tmu || is2te) && vlbe.size()==2 && vtaue.size()==2 ){

vj0.SetPxPyPzE(vjpx[0], vjpy[0], vjpz[0], vje[0]);
vj1.SetPxPyPzE(vjpx[1], vjpy[1], vjpz[1], vje[1]);
//use svfit to have correct mass

//cat 1
if(ngj==2){

if( (vj0+vj1).M()<140 && (vj0+vj1).M()>40 )  //w-tag
c1sl++;
//cat 2
else 
c2sl++;
}
else if(ngj==1){

vb0.SetPxPyPzE(vlbpx[0], vlbpy[0], vlbpz[0], vlbe[0]);
//cat 3
if( (65<(vb0+vj0).M()<95) || (65<(vb0+vj1).M()<95) )
c3sl++;
//cat 4
else 
c4sl++;
}

//cat 5
else if(ngj>=3){
c5sl++;
}
}

//////////////////////////
// event categories: DL //
//////////////////////////

if( (is2tmu || is2te || istmue || istlmu || istle || istmule || istelmu) && vtaue.size()==2  ){

//cat 6
if( vtbe.size()==2 ){
c6dl++;
}
//cat 7
else if( vtbe.size()==2 && vlbe.size()==1  ){
c7dl++;
}

}
 */

}//evt


//if( (sID&1<<0)==1<<0 ) weight = (20.*1000.)/*lumi 20fb-1*/*(0.1293/*tth*/*0.0632/*htautau*/)/999997.; //tth  
//if( (sID&1<<2)==1<<2 ) weight = (20.*1000.)*234.7/4908667.; //ttjets  

cout << "Nentries: " <<nentries << ", run on: "<<nev <<" events =? " << nmuORe << endl;

}
