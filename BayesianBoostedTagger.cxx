// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include "BayesianBoostedTagger.h"
#include <BAT/BCSummaryTool.h>

#include "TF1.h"

#include <BAT/BCMath.h>

// ---------------------------------------------------------
BayesianBoostedTagger::BayesianBoostedTagger() : BCModel(), 
						 m_runMode(kOmicron),
						 _doDiagnostics(false),
						 _Nevent(0), _max_E(1*TeV), _max_eta(5.0), 
						 m_calo(0), m_template(0), 
						 _t_truth(NULL), _b_truth(NULL), _W_truth(NULL), _q1_truth(NULL), _q2_truth(NULL)
{
   // default constructor
   DefineParameters();
}

// ---------------------------------------------------------
BayesianBoostedTagger::BayesianBoostedTagger(const char * name, BBTRunMode mode ) : BCModel(name), 
										    m_runMode(mode),
										    _doDiagnostics(false),
										    _Nevent(0),_max_E(1*TeV), _max_eta(2.5), 
										    m_calo(0), m_template(0), 
										    _t_truth(NULL), _b_truth(NULL), _W_truth(NULL), _q1_truth(NULL), _q2_truth(NULL)
{
   // constructor
   DefineParameters();
}

// ---------------------------------------------------------
BayesianBoostedTagger::~BayesianBoostedTagger()
   // default destructor
{
  if( m_calo ) delete m_calo;
  if( m_template ) delete m_template;
  // if( m_pTfunc ) delete m_pTfunc;
}

// ---------------------------------------------------------
void BayesianBoostedTagger::SetJet( const fastjet::ClusterSequence & cs,
				    const fastjet::PseudoJet & jet )
{
  _cs  = &cs;
  _jet = fastjet::PseudoJet( jet );
}


// ---------------------------------------------------------
void BayesianBoostedTagger::DefineParameters()
{
   // Add parameters to your model here.
   // You can then use them in the methods below by calling the
   // parameters.at(i) or parameters[i], where i is the index
   // of the parameter. The indices increase from 0 according to the
   // order of adding the parameters.

  AddParameter( "t_pt", 0., 1.5*TeV );
  if( m_runMode == kOmega ) {
    AddParameter( "t_eta", -2.5, 2.5 );
    AddParameter( "t_phi", -3.14, 3.14 );
  }

  /*
  AddParameter( "tWb_theta", 0, M_PI );
  AddParameter( "tWb_phi", -M_PI, M_PI );

  AddParameter( "Wqq_theta", 0, M_PI );
  AddParameter( "Wqq_phi", -M_PI, M_PI );
  */

  AddParameter( "tWb_theta", 0, 1 );
  AddParameter( "tWb_phi",   0, 1 );
  AddParameter( "Wqq_theta", 0,	1 );
  AddParameter( "Wqq_phi",   0, 1 );

  //TF1 * pTfunc = new TF1( "pTfunc", "exp( -x/100. )", 0, _max_E );
//  TF1 * pTfunc = new TF1( "pTfunc", "1 - x/1000 ", 0, _max_E );
//  TF1 * pTfunc = new TF1( "pTfunc", "TMath::Landau(x,[0],[1],0)", 0, _max_E );
//  pTfunc->SetParameters( 50, _max_E/10. );
  // m_pTfunc->SetParName(0, "tau" );

  //m_pTfunc->SetParameter( 0, _max_E );

  SetPriorConstantAll();

//  delete pTfunc;

//     MCMCSetPrecision( BCEngineMCMC::kLow );
//    MCMCSetPrecision( BCEngineMCMC::kMedium );
//  MCMCSetPrecision( BCEngineMCMC::kHigh );

/*  
  MCMCSetNChains(10); 
  MCMCSetNIterationsRun(20000); 
  MCMCSetNIterationsMax(10000); 
  MCMCSetNIterationsUpdate(1000);
 */

  SetIntegrationMethod( BCModel::kIntMonteCarlo );
  //SetIntegrationMethod( BCModel::kIntMetropolis );
  //SetIntegrationMethod( BCModel::kIntImportance );

  // 5 DoF
  if( m_runMode == kOmicron ) {
    MCMCSetNChains(5); 
    MCMCSetNIterationsRun(100000); 
    MCMCSetNIterationsMax(100000); 
    MCMCSetNIterationsUpdate(10000);
  }
  // 7 DoF
  else if( m_runMode == kOmega ) {
    MCMCSetNChains(5); 
    MCMCSetNIterationsRun(500000); 
    MCMCSetNIterationsMax(500000); 
    MCMCSetNIterationsUpdate(5000);
  }

  

  const int n_eta_cells = 200;
  const int n_phi_cells = 128;
  
//  const int n_eta_cells = 50;
//  const int n_phi_cells = 64;

  m_calo     = new TH2D( "calo",     "Calo",     n_eta_cells, -_max_eta, _max_eta, n_phi_cells, -M_PI, M_PI );
  m_template = new TH2D( "template", "Template", n_eta_cells, -_max_eta, _max_eta, n_phi_cells, -M_PI, M_PI );

  const int    lhood_nbins = 200;
  const double lhood_max   = 10;
  const int    eta_nbins   = 100;
  const int    phi_nbins   = 64;
  const int    E_nbins     = 150;
  const double E_max       = 300;

  m_lhood = new TH1D( "lhood", "-logL", lhood_nbins, 0., lhood_max );
  
  m_t_pT_vs_lhood = new TH2D( "t_pT_vs_lhood", "p_{T}(t) vs -logL", lhood_nbins, 0., lhood_max, 100, 0., 2*TeV );

  m_h_theta_bW_vs_lhood = new TH2D( "theta_bW_vs_lhood", "#theta(b,W) vs lhood",   lhood_nbins, 0, lhood_max, 100, 0, M_PI );
  m_h_theta_qq_vs_lhood = new TH2D( "theta_qq_vs_lhood", "#theta(q1,q2) vs lhood", lhood_nbins, 0, lhood_max, 100, 0, M_PI );

  m_h_W_pT_vs_lhood    = new TH2D( "W_pT_vs_lhood",  "W pT vs lhood",   lhood_nbins, 0, lhood_max, E_nbins, 0., E_max );
  m_h_W_eta_vs_lhood   = new TH2D( "W_eta_vs_lhood", "W #eta vs lhood", lhood_nbins, 0, lhood_max, eta_nbins, -2.5, 2.5 );
  m_h_W_phi_vs_lhood   = new TH2D( "W_phi_vs_lhood", "W #phi vs lhood", lhood_nbins, 0, lhood_max, phi_nbins, -M_PI, M_PI );
  m_h_W_E_vs_lhood     = new TH2D( "W_E_vs_lhood",   "W E vs lhood",    lhood_nbins, 0, lhood_max, E_nbins, 0., E_max );

  m_h_b_pT_vs_lhood    = new TH2D( "b_pT_vs_lhood",  "b pT vs lhood",   lhood_nbins, 0, lhood_max, E_nbins, 0., E_max );
  m_h_b_eta_vs_lhood   = new TH2D( "b_eta_vs_lhood", "b #eta vs lhood", lhood_nbins, 0, lhood_max, eta_nbins, -2.5, 2.5 );
  m_h_b_phi_vs_lhood   = new TH2D( "b_phi_vs_lhood", "b #phi vs lhood", lhood_nbins, 0, lhood_max, phi_nbins, -M_PI, M_PI );
  m_h_b_E_vs_lhood     = new TH2D( "b_E_vs_lhood",   "b E vs lhood",    lhood_nbins, 0, lhood_max, E_nbins, 0., E_max );

  m_h_q1_pT_vs_lhood    = new TH2D( "q1_pT_vs_lhood",  "q1 pT vs lhood",   lhood_nbins, 0, lhood_max, E_nbins, 0., E_max );
  m_h_q1_eta_vs_lhood   = new TH2D( "q1_eta_vs_lhood", "q1 #eta vs lhood", lhood_nbins, 0, lhood_max, eta_nbins, -2.5, 2.5 );
  m_h_q1_phi_vs_lhood   = new TH2D( "q1_phi_vs_lhood", "q1 #phi vs lhood", lhood_nbins, 0, lhood_max, phi_nbins, -M_PI, M_PI );
  m_h_q1_E_vs_lhood     = new TH2D( "q1_E_vs_lhood",   "q1 E vs lhood",    lhood_nbins, 0, lhood_max, E_nbins, 0., E_max );

  m_h_q2_pT_vs_lhood    = new TH2D( "q2_pT_vs_lhood",  "q2 pT vs lhood",   lhood_nbins, 0, lhood_max, E_nbins, 0., E_max );
  m_h_q2_eta_vs_lhood   = new TH2D( "q2_eta_vs_lhood", "q2 #eta vs lhood", lhood_nbins, 0, lhood_max, eta_nbins, -2.5, 2.5 );
  m_h_q2_phi_vs_lhood   = new TH2D( "q2_phi_vs_lhood", "q2 #phi vs lhood", lhood_nbins, 0, lhood_max, phi_nbins, -M_PI, M_PI );
  m_h_q2_E_vs_lhood     = new TH2D( "q2_E_vs_lhood",   "q2 E vs lhood",    lhood_nbins, 0, lhood_max, E_nbins, 0., E_max );

  m_h_Rqq_vs_lhood      = new TH2D( "Rqq_vs_lhood",    "#Delta R(q,q) vs lhood", lhood_nbins, 0, lhood_max, 100, 0, 5 );
  m_h_RbW_vs_lhood      = new TH2D( "RbW_vs_lhood",    "#Delta R(b,W) vs lhood", lhood_nbins, 0, lhood_max, 100, 0, 5 );

  m_h_t_truth_dR_vs_lhood   = new TH2D( "t_truth_dR_vs_lhood",   "#Delta R(t_{reco}, t_{truth})",    lhood_nbins, 0, lhood_max, 100, 0, 10 );
  m_h_t_truth_dPt_vs_lhood  = new TH2D( "t_truth_dPt_vs_lhood",  "#Delta pT(t_{reco}, t_{truth})",   lhood_nbins, 0, lhood_max, 100, -1, 1 );
  m_h_t_truth_dEta_vs_lhood = new TH2D( "t_truth_dEta_vs_lhood", "#Delta #eta(t_{reco}, t_{truth})", lhood_nbins, 0, lhood_max, 100, -1, 1 );
  m_h_t_truth_dPhi_vs_lhood = new TH2D( "t_truth_dPhi_vs_lhood", "#Delta #phi(t_{reco}, t_{truth})", lhood_nbins, 0, lhood_max, 100, -1, 1 );

  m_h_b_truth_dR_vs_lhood   = new TH2D( "b_truth_dR_vs_lhood",   "#Delta R(b_{reco}, b_{truth})",    lhood_nbins, 0, lhood_max, 100, 0, 10 );
  m_h_b_truth_dPt_vs_lhood  = new TH2D( "b_truth_dPt_vs_lhood",  "#Delta pT(b_{reco}, b_{truth})",   lhood_nbins, 0, lhood_max, 100, -1, 1 );
  m_h_b_truth_dEta_vs_lhood = new TH2D( "b_truth_dEta_vs_lhood", "#Delta #eta(b_{reco}, b_{truth})", lhood_nbins, 0, lhood_max, 100, -1, 1 );
  m_h_b_truth_dPhi_vs_lhood = new TH2D( "b_truth_dPhi_vs_lhood", "#Delta #phi(b_{reco}, b_{truth})", lhood_nbins, 0, lhood_max, 100, -1, 1 );

  m_h_W_truth_dR_vs_lhood   = new TH2D( "W_truth_dR_vs_lhood",   "#Delta R(W_{reco}, W_{truth})",    lhood_nbins, 0, lhood_max, 100, 0, 10 );
  m_h_W_truth_dPt_vs_lhood  = new TH2D( "W_truth_dPt_vs_lhood",  "#Delta pT(W_{reco}, W_{truth})",   lhood_nbins, 0, lhood_max, 100, -1, 1 );
  m_h_W_truth_dEta_vs_lhood = new TH2D( "W_truth_dEta_vs_lhood", "#Delta #eta(W_{reco}, W_{truth})", lhood_nbins, 0, lhood_max, 100, -1, 1 );
  m_h_W_truth_dPhi_vs_lhood = new TH2D( "W_truth_dPhi_vs_lhood", "#Delta #phi(W_{reco}, W_{truth})", lhood_nbins, 0, lhood_max, 100, -1, 1 );
}

#define PRINTOUT(x) printf( "Event %i / ( %4.1f, %3.2f, %3.2f, %4.1f ; %3.1f )\n", _Nevent, x.Pt(), x.Eta(), x.Phi(), x.E(), x.M() );
// ---------------------------------------------------------
void BayesianBoostedTagger::MakeParticles( const vector<double>& parameters )
{
  const double t_pt  = parameters.at(0);  
  const double t_eta = ( m_runMode == kOmega ) ? parameters.at(1) : _fjet.Eta();
  const double t_phi = ( m_runMode == kOmega ) ? parameters.at(2) : _fjet.Phi();
  _t.SetPtEtaPhiM( t_pt, t_eta, t_phi, masst );

  //PRINTOUT(_t);

  //////////////////
  // decay t->Wb

  double W_E = ( masst*masst + massW*massW - massb*massb ) / ( 2. * masst );
  double b_E = ( masst*masst - massW*massW + massb*massb ) / ( 2. * masst );
  
  double p = sqrt( W_E*W_E - massW*massW ); // conservation of momentum p_W = p_b

  // isotropic decay angles
  int idx = ( m_runMode == kOmicron ) ? 1 : 3;
  const double tWb_theta = acos( 2 * parameters.at(idx) - 1 ); 
  idx = ( m_runMode == kOmicron ) ? 2 : 4;
  const double tWb_phi   = 2 * M_PI * parameters.at(idx);

  double px = p * sin(tWb_theta) * cos(tWb_phi);
  double py = p * sin(tWb_theta) * sin(tWb_phi);
  double pz = p * cos(tWb_theta);
  _W.SetPxPyPzE( px, py, pz, W_E );
  _b.SetPxPyPzE( -px, -py, -pz, b_E );
  
  // boost
  Boost( _W, _t );
  Boost( _b, _t );
  // _W.Boost( _t.Vect() );
  //_b.Boost( _t.Vect() );

  //PRINTOUT( _W );
  //PRINTOUT( _b );

  //////////////////
  // decay W->qq
  // assume both q massless

  double q1_E = massW / 2;
  double q2_E = massW / 2;
  double Wqq_P = massW / 2.;
    

  idx = ( m_runMode == kOmicron ) ? 3 : 5;
  const double Wqq_theta = acos( 2 * parameters.at(idx) - 1 );
  idx = ( m_runMode == kOmicron ) ? 4 : 6;
  const double Wqq_phi   = 2 * M_PI * parameters.at(idx);


  px = Wqq_P * sin(Wqq_theta) * cos(Wqq_phi);
  py = Wqq_P * sin(Wqq_theta) * sin(Wqq_phi);
  pz = Wqq_P * cos(Wqq_theta);
  _q1.SetPxPyPzE( px, py, pz, q1_E );
  _q2.SetPxPyPzE( -px, -py, -pz, q2_E );
  // _q1.SetPxPyPzM( px, py, pz, 0. );
  //_q2.SetPxPyPzM( -px, -py, -pz, 0. );
  
  // boost
  Boost( _q1, _W );
  Boost( _q2, _W );
  //_q1.Boost( _W.Vect() );
  //_q2.Boost( _W.Vect() );

  //PRINTOUT( _q1 );
  // PRINTOUT( _q2 );
}


// ---------------------------------------------------------
void BayesianBoostedTagger::Boost( TLorentzVector& p, const TLorentzVector& parent )
{
  // beta and gamma values
  double betax = parent.Px() / parent.E();
  double betay = parent.Py() / parent.E();
  double betaz = parent.Pz() / parent.E();
  double beta2 = betax*betax + betay*betay + betaz*betaz;
  double gamma = 1.0/sqrt(1.0-beta2);
  double dot   = betax*p.Px() + betay*p.Py() + betaz*p.Pz();
  double prod  = gamma*( gamma*dot/(1.0+gamma) + p.E() );

  double pX = p.Px() + betax*prod;
  double pY = p.Py() + betay*prod;
  double pZ = p.Pz() + betaz*prod;
  double e  = gamma*(p.E() + dot);

  p.SetPxPyPzE( pX, pY, pZ, e );
}


// ---------------------------------------------------------
void BayesianBoostedTagger::SetTruthPartons( const TLorentzVector& t, const TLorentzVector& b, const TLorentzVector& W, const TLorentzVector& q1, const TLorentzVector& q2 )
{
  _t_truth  = &t;
  _b_truth  = &b;
  _W_truth  = &W;
  _q1_truth = &q1;
  _q2_truth = &q2;
}


// ---------------------------------------------------------
double BayesianBoostedTagger::RunTagger( const TH2D& hcalo )
{
   m_calo->Reset();
   m_calo->Add( &hcalo );
   return RunTagger();
}

// ---------------------------------------------------------
double BayesianBoostedTagger::RunTagger( const vector<TLorentzVector>& clusters, const TLorentzVector * fjet, bool doDiagnostics )
{
  _doDiagnostics = doDiagnostics;

   const size_t n_cl = clusters.size();

   m_clusters.clear();
   // m_clusters.reserve( n_cl );
   // copy( clusters.begin(), clusters.end(), m_clusters.begin() );
   m_clusters = clusters;

   for( size_t k = 0 ; k < n_cl ; ++k ) {
     const double eta = clusters.at(k).Eta();
     const double phi = clusters.at(k).Phi();
     const double E   = clusters.at(k).E();

     // cout << "cluster " << k << " " << eta << " " << phi << " " << E << endl;
     m_calo->Fill( eta, phi, E );
   }
   
   if( fjet ) {
     _fjet = *fjet;
   }
   else {
     cout << "WARNING: no fat jet defined" << endl;
   }

   const double lhood = RunTagger(); 

   return lhood;
}


// ---------------------------------------------------------
double BayesianBoostedTagger::RunTagger( const fastjet::PseudoJet & jet )
{
  _jet = jet;

  TLorentzVector fjet;
  fjet.SetPtEtaPhiE( jet.perp(), jet.eta(), jet.phi() - 3.14, jet.E() );

  vector< TLorentzVector > clusters;
  const vector<fastjet::PseudoJet> jc = _jet.constituents();
  
  for( size_t k = 0 ; k < jc.size() ; ++k ) {
    const double pt  = jc.at(k).perp();
    const double eta = jc.at(k).eta();
    const double phi = jc.at(k).phi() - 3.14;
    const double E   = jc.at(k).E();

    TLorentzVector cl;
    cl.SetPtEtaPhiE( pt, eta, phi, E );
    clusters.push_back( cl );

    // cout << "cluster " << k << " " << eta << " " << phi << " " << E << endl;
    m_calo->Fill( eta, phi, E );
  }

  return RunTagger( clusters, &fjet );
}

// ---------------------------------------------------------
double BayesianBoostedTagger::RunTagger()
{
  // first reset diagnostics
  m_lhood->Reset();

  m_t_pT_vs_lhood->Reset();
 
  m_h_W_pT_vs_lhood->Reset();
  m_h_W_eta_vs_lhood->Reset();
  m_h_W_phi_vs_lhood->Reset();
  m_h_W_E_vs_lhood->Reset();

  m_h_b_pT_vs_lhood->Reset();
  m_h_b_eta_vs_lhood->Reset();
  m_h_b_phi_vs_lhood->Reset();
  m_h_b_E_vs_lhood->Reset();

  m_h_q1_pT_vs_lhood->Reset();
  m_h_q1_eta_vs_lhood->Reset();
  m_h_q1_phi_vs_lhood->Reset();
  m_h_q1_E_vs_lhood->Reset();

  m_h_q2_pT_vs_lhood->Reset();
  m_h_q2_eta_vs_lhood->Reset();
  m_h_q2_phi_vs_lhood->Reset();
  m_h_q2_E_vs_lhood->Reset();

  m_h_Rqq_vs_lhood->Reset();
  m_h_RbW_vs_lhood->Reset();

  m_h_t_truth_dR_vs_lhood->Reset();
  m_h_t_truth_dPt_vs_lhood->Reset();
  m_h_t_truth_dEta_vs_lhood->Reset();
  m_h_t_truth_dPhi_vs_lhood->Reset();
  
  m_h_b_truth_dR_vs_lhood->Reset();
  m_h_b_truth_dPt_vs_lhood->Reset();
  m_h_b_truth_dEta_vs_lhood->Reset();
  m_h_b_truth_dPhi_vs_lhood->Reset();

  m_h_W_truth_dR_vs_lhood->Reset();
  m_h_W_truth_dPt_vs_lhood->Reset();
  m_h_W_truth_dEta_vs_lhood->Reset();
  m_h_W_truth_dPhi_vs_lhood->Reset();
  

  // if( m_calo->GetEntries() == 0 ) FillCalo();

  /*
  //Normalize();
  MarginalizeAll();
  FindMode(); // GetBestFitParameters() );

  //  const vector<double> bfp = GetBestFitParametersMarginalized();
  const vector<double> bfp = GetBestFitParameters();
  // cout << bfp[0] << " " << bfp[1] << " " << bfp[2] << " " << bfp[3] << " " << bfp[4] << endl;
  //if( bfp.size() == 0 ) cout << "ERROR: no best fit parameters" << endl;
  
  //FindMode( bfp );  
  //FindModeMinuit( bfp, -1 );

  //bfp = GetBestFitParameters();
  //cout << bfp[0] << " " << bfp[1] << " " << bfp[2] << " " << bfp[3] << " " << bfp[4] << endl;

  MakeParticles( bfp );

  const double lhood = LogLikelihood( bfp );
  //const double lhood = LogProbability( bfp ); 
  //const double lhood = log( GetNormalization() );
  */
  
  MarginalizeAll();
  const vector<double> bfp = GetBestFitParameters();
  //FindModeMinuit( bfp, -1 );
  FindMode( bfp );
  const int status = GetMinuitErrorFlag();
  printf("DEBUG: minuit status %i\n", status );
  const double lhood = LogLikelihood( bfp );

  ++_Nevent;

  return lhood;
}

// ---------------------------------------------------------
void BayesianBoostedTagger:: Dump( const char *	tag )
{
  BCSummaryTool * summary = new BCSummaryTool( this );
  char buf[128];
  //sprintf( buf, "bbt_%s_update.ps", tag );
  //summary->PrintKnowledgeUpdatePlots( buf );

  sprintf( buf, "bbt_%s_plots.ps", tag );
  PrintAllMarginalized( buf );

  //sprintf( buf, "bbt_%i_correlation.ps", _Nevent );
  //summary->PrintCorrelationPlot( buf );

  TFile * fdiagnostics = TFile::Open( "bbt_diagnostics.root", "update" );
  fdiagnostics->cd();
  fdiagnostics->mkdir( tag );
  fdiagnostics->cd( tag );

  m_lhood->Write();

  m_t_pT_vs_lhood->Write();

  m_h_W_pT_vs_lhood->Write();
  m_h_W_eta_vs_lhood->Write();
  m_h_W_phi_vs_lhood->Write();
  m_h_W_E_vs_lhood->Write();

  m_h_b_pT_vs_lhood->Write();
  m_h_b_eta_vs_lhood->Write();
  m_h_b_phi_vs_lhood->Write();
  m_h_b_E_vs_lhood->Write();

  m_h_q1_pT_vs_lhood->Write();
  m_h_q1_eta_vs_lhood->Write();
  m_h_q1_phi_vs_lhood->Write();
  m_h_q1_E_vs_lhood->Write();

  m_h_q2_pT_vs_lhood->Write();
  m_h_q2_eta_vs_lhood->Write();
  m_h_q2_phi_vs_lhood->Write();
  m_h_q2_E_vs_lhood->Write();

  m_h_Rqq_vs_lhood->Write();
  m_h_RbW_vs_lhood->Write();

  m_calo->Write();
  m_template->Write();

  m_h_theta_bW_vs_lhood->Write();
  m_h_theta_qq_vs_lhood->Write();

  m_h_t_truth_dR_vs_lhood->Write();
  m_h_t_truth_dPt_vs_lhood->Write();
  m_h_t_truth_dEta_vs_lhood->Write();
  m_h_t_truth_dPhi_vs_lhood->Write();

  m_h_b_truth_dR_vs_lhood->Write();
  m_h_b_truth_dPt_vs_lhood->Write();
  m_h_b_truth_dEta_vs_lhood->Write();
  m_h_b_truth_dPhi_vs_lhood->Write();

  m_h_W_truth_dR_vs_lhood->Write();
  m_h_W_truth_dPt_vs_lhood->Write();
  m_h_W_truth_dEta_vs_lhood->Write();
  m_h_W_truth_dPhi_vs_lhood->Write();

  fdiagnostics->cd();
  fdiagnostics->Close();

  delete fdiagnostics;
  delete summary;
}

// ---------------------------------------------------------
void BayesianBoostedTagger::FillCalo()
{
  m_calo->Reset();

  const vector<fastjet::PseudoJet> jc = _jet.constituents();
  
  for( size_t k = 0 ; k < jc.size() ; ++k ) {
    const double eta = jc.at(k).eta();
    const double phi = jc.at(k).phi() - 3.14;
    const double E   = jc.at(k).E();

    // cout << "cluster " << k << " " << eta << " " << phi << " " << E << endl;
    m_calo->Fill( eta, phi, E );
  }
}

// ---------------------------------------------------------
double BayesianBoostedTagger::LogLikelihood(const std::vector<double> &parameters)
{
   // This methods returns the logarithm of the conditional probability
   // p(data|parameters). This is where you have to define your model.

   double logprob = 0.;

   MakeParticles( parameters );
 
   vector<TLorentzVector> pv;
   pv.push_back( _b );
   pv.push_back( _q1 );
   pv.push_back( _q2 );

   logprob += LogOverlap( pv ); 

   return logprob;
}


// ---------------------------------------------------------

/*
double BayesianBoostedTagger::LogAPrioriProbability(const std::vector<double> &parameters)
{
   // This method returns the logarithm of the prior probability for the
   // parameters p(parameters).

   double logprob = 0.;
   //  return logprob;
//   double x = parameters.at(0);
//   double y = parameters.at(1);

   return logprob;
}
*/
// ---------------------------------------------------------

// ---------------------------------------------------------
void BayesianBoostedTagger::MCMCIterationInterface()
{
  if( !_doDiagnostics ) return;

  int nchains = this -> MCMCGetNChains(); 
  
  for( int i = 0; i < nchains; ++i ) {
    vector< double > parameters = MCMCGetx(i);

    const double logL = LogLikelihood( parameters );
    //const double logL = LogProbability( parameters );
    // const double logL = this -> MCMCGetLogProbx(i);
    
    const double t_pt = _t.Pt();

    m_lhood->Fill( -logL );

    m_t_pT_vs_lhood->Fill( -logL, t_pt );

    m_h_W_pT_vs_lhood->Fill(  -logL, _W.Pt()  );
    m_h_W_eta_vs_lhood->Fill( -logL, _W.Eta() );
    m_h_W_phi_vs_lhood->Fill( -logL, _W.Phi() );
    m_h_W_E_vs_lhood->Fill(   -logL, _W.E()   );

    m_h_b_pT_vs_lhood->Fill(  -logL, _b.Pt()  );
    m_h_b_eta_vs_lhood->Fill( -logL, _b.Eta() );
    m_h_b_phi_vs_lhood->Fill( -logL, _b.Phi() );
    m_h_b_E_vs_lhood->Fill(   -logL, _b.E()   );

    m_h_q1_pT_vs_lhood->Fill(  -logL, _q1.Pt()  );
    m_h_q1_eta_vs_lhood->Fill( -logL, _q1.Eta() );
    m_h_q1_phi_vs_lhood->Fill( -logL, _q1.Phi() );
    m_h_q1_E_vs_lhood->Fill(   -logL, _q1.E()   );

    m_h_q2_pT_vs_lhood->Fill(  -logL, _q2.Pt()  );
    m_h_q2_eta_vs_lhood->Fill( -logL, _q2.Eta() );
    m_h_q2_phi_vs_lhood->Fill( -logL, _q2.Phi() );
    m_h_q2_E_vs_lhood->Fill(   -logL, _q2.E()   );

    const double theta_bW = _W.Angle( _b.Vect() );
    const double theta_qq = _q1.Angle( _q2.Vect() );

    m_h_theta_bW_vs_lhood->Fill( -logL, theta_bW );
    m_h_theta_qq_vs_lhood->Fill( -logL, theta_qq );

    const double dRqq = _q1.DeltaR( _q2 );
    m_h_Rqq_vs_lhood->Fill( -logL, dRqq );
    
    const double dRbW = _b.DeltaR( _W );
    m_h_RbW_vs_lhood->Fill( -logL, dRbW );

    // now compare to truth, if available
    if( _t_truth ) {
      const double dR   = _t_truth->DeltaR( _t );
      const double dPt  = ( _t.Pt() - _t_truth->Pt() )  / _t_truth->Pt();
      const double dEta = ( _t.Eta() - _t_truth->Eta() ) / _t_truth->Eta();
      const double dPhi = ( _t_truth->DeltaPhi( _t ) ) / _t_truth->Phi() ;

      m_h_t_truth_dR_vs_lhood->Fill( -logL, dR );
      m_h_t_truth_dPt_vs_lhood->Fill( -logL, dPt );
      m_h_t_truth_dEta_vs_lhood->Fill( -logL, dEta );
      m_h_t_truth_dPhi_vs_lhood->Fill( -logL, dPhi );
    }
    if( _b_truth ) {
      const double dR   = _b_truth->DeltaR( _b );
      const double dPt  = ( _b.Pt() - _b_truth->Pt() )  / _b_truth->Pt();
      const double dEta = ( _b.Eta() - _b_truth->Eta() ) / _b_truth->Eta();
      const double dPhi = ( _b_truth->DeltaPhi( _b ) ) / _b_truth->Phi() ;

      m_h_b_truth_dR_vs_lhood->Fill( -logL, dR );
      m_h_b_truth_dPt_vs_lhood->Fill( -logL, dPt );
      m_h_b_truth_dEta_vs_lhood->Fill( -logL, dEta );
      m_h_b_truth_dPhi_vs_lhood->Fill( -logL, dPhi );
    }
    if( _W_truth ) {
      const double dR   = _W_truth->DeltaR( _W );
      const double dPt  = ( _W.Pt() - _W_truth->Pt() )  / _W_truth->Pt();
      const double dEta = ( _W.Eta() - _W_truth->Eta() ) / _W_truth->Eta();
      const double dPhi = ( _W_truth->DeltaPhi( _W ) ) / _W_truth->Phi() ;

      m_h_W_truth_dR_vs_lhood->Fill( -logL, dR );
      m_h_W_truth_dPt_vs_lhood->Fill( -logL, dPt );
      m_h_W_truth_dEta_vs_lhood->Fill( -logL, dEta );
      m_h_W_truth_dPhi_vs_lhood->Fill( -logL, dPhi );
    }
  }

}


// ---------------------------------------------------------
double BayesianBoostedTagger::LogOverlap( const vector<TLorentzVector>& pv ) const
{
  double ov = 0;

  const size_t np = pv.size();
/*
  m_template->Reset();
  for( size_t a = 0 ; a < pv.size() ; ++a ) {
    const TLorentzVector * p = &pv.at(a);

    const double E_a = p->E();
    const double eta_a = p->Eta();
    const double phi_a = p->Phi();

    m_template->Fill( eta_a, phi_a, E_a );
  }
*/
  for( size_t a = 0 ; a < pv.size() ; ++a ) {
    const TLorentzVector * p = &pv.at(a);

    //    const double eta_a = p->Eta();
    // const double phi_a = p->Phi();

    //  const int bin_eta = m_template->GetXaxis()->FindBin( eta_a );
    //  const int bin_phi = m_template->GetYaxis()->FindBin( phi_a );

//    const double E_p = m_template->GetBinContent( bin_eta, bin_phi );
    const double E_p = p->E();
    //           const double sigma_p = E_p / 2.; // template tagging paper
        const double sigma_p = E_p / 3.;
    //    const double sigma_p = 0.1 * E_p; // avg ATLAS
//    const double sigma_p = 20.;
//    const double sigma_p = E_p * 0.2 * exp( -E_p / 500. ); 
//    const double sigma_p = 10. + 0.0005*E_p*E_p;
//    const double sigma_p =  E_p * ( 0.05 + 0.5/sqrt(E_p) + 0.1/E_p ); // fit ATLAS Jet Energy Resolution
//    const double sigma_p =  E_p * ( 1.4/sqrt(E_p) - 2.7/E_p );
//    const double sigma_p = 15.;

    double E_calo = 0.; 

//    E_calo = m_calo->Integral( bin_eta-1, bin_eta+1, bin_phi-1, bin_phi+1 );
    //E_calo = m_calo->GetBinContent( bin_eta, bin_phi );
    
    if( m_clusters.size() == 0 ) {
      cout << "ERROR: no clusters in fat jet?" << endl;
    }

/*
    // dR =  2M / pT
    double dR_max = p->M() / ( 2. * p->E() );
    dR_max = ( dR_max > 0.3 ) ? 0.3 : dR_max;
    dR_max = ( dR_max < 0.1 ) ? 0.1 : dR_max;
*/
   double dR_max = 0.2;
    //   double dR_max = ( p->Pt() < 40. ) ? 0.35 : 100 / p->Pt();
    //cout << "parton pT = " << p->Pt() << " m = " << p->M() << " dR_max = " << dR_max << endl;

    TLorentzVector calo;
    for( vector<TLorentzVector>::const_iterator i_cl = m_clusters.begin() ; i_cl != m_clusters.end() ; ++i_cl ) {
      TLorentzVector cl( (*i_cl) );

       const double dR = p->DeltaR( cl );
       if( dR > dR_max ) continue;
  
       E_calo += cl.E();
       //    calo += cl;
    } 

    // if( dE != 0. ) cout << "parton " << a << " dE = " << dE << " E_p = " << E_p << endl;

    // LogGauss 
//    const double dE = E_calo - E_p;


   if( E_calo < 10. ) {
      //ov += -HUGELHOOD;
      //continue;
      return -HUGELHOOD;
    }

   const double alpha = 1.0;
   const double dE = ( alpha * E_p - E_calo ) / sigma_p;
   ov -= dE * dE;

    // cout << "E_p = " << E_p << " E_calo = " << calo.E() << " ov = " << ov << endl;

    //    E_calo = calo.E();
   // ov += BCMath::LogGaus( E_calo, E_p, sigma_p, true );
   //  ov += BCMath::LogGaus( E_calo, E_p, sigma_p, false );
  }
  // cout << ov << endl;

  return ov;
}


void BayesianBoostedTagger::Compare2Truth()
{
  /*
  const double dR   = _t_truth.DeltaR( _t );
  const double dPhi = _t_truth.DeltaPhi( _t );
  const double dPt  = _t_truth.Pt() - _t.Pt();
  */
}
