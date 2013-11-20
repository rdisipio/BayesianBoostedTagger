// ***************************************************************
// This file was created using the CreateProject.sh script.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#ifndef __BAT__BAYESIANBOOSTEDTAGGER__H
#define __BAT__BAYESIANBOOSTEDTAGGER__H

#include <iterator>
#include <vector>
#include <math.h>

#include <BAT/BCModel.h>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>

#include "TFile.h"
#include "TLorentzVector.h"
#include "TH2D.h"
#include "TF1.h"

#define HUGELHOOD 1e20
#define TeV 1000.
#define massW  80.4
#define masst  172.5
#define massb  4.5
#define gammaW 2.1
#define gammaTop 1.5

// This is a BayesianBoostedTagger header file.
// Model source code is located in file BayesianBoostedTagger/BayesianBoostedTagger.cxx

using namespace std;

// ---------------------------------------------------------
class BayesianBoostedTagger : public BCModel
{
   public:
      enum BBTRunMode { kOmicron, kOmega }; 

   public:

      // Constructors and destructor
      BayesianBoostedTagger();
      BayesianBoostedTagger(const char * name, BBTRunMode mode = kOmicron );
      ~BayesianBoostedTagger();

      void SetJet( const fastjet::ClusterSequence & cs,
		 const fastjet::PseudoJet & jet );

      void SetTruthPartons( const TLorentzVector& t, const TLorentzVector& b, const TLorentzVector& W, const TLorentzVector& q1, const TLorentzVector& q2 );  

      double RunTagger( const fastjet::PseudoJet & jet );
      double RunTagger( const TH2D& hcalo );
      double RunTagger( const vector<TLorentzVector>& clusters, const TLorentzVector * fjet = NULL, bool doDiagnostics = false );
      double RunTagger();
      void   Dump( const char * tag = "bbt" );

      void FillCalo();

      // Methods to overload, see file BayesianBoostedTagger.cxx
      void DefineParameters();
      //double LogAPrioriProbability(const std::vector<double> &parameters);
      double LogLikelihood(const std::vector<double> &parameters);
      void MCMCIterationInterface();

      void   MakeParticles( const vector<double>& parameters );
      double LogOverlap( const vector< TLorentzVector >& pv ) const;

      inline TLorentzVector GetTopCandidate() { return _t; };
      inline TLorentzVector GetWCandidate()   { return _W; };
      inline TLorentzVector GetBCandidate()   { return _b; };
      inline TLorentzVector GetQ1Candidate()  { return _q1; };
      inline TLorentzVector GetQ2Candidate()  { return _q2; };

      void Boost( TLorentzVector& p,  const TLorentzVector& parent );
      
      inline void SetRunMode( BBTRunMode mode ) { m_runMode = mode; };

      void Compare2Truth(); 

 private:
      BBTRunMode m_runMode;
      bool  _doDiagnostics;

      int    _Nevent;
      double _max_E;
      double _max_eta;

//      TF1  * m_pTfunc;
      vector< TLorentzVector > m_clusters;
      TH2D * m_calo;
      TH2D * m_template;
      TH1D * m_lhood;

      TH2D * m_t_pT_vs_lhood;

      TH2D * m_h_W_pT_vs_lhood;
      TH2D * m_h_W_eta_vs_lhood;
      TH2D * m_h_W_phi_vs_lhood;
      TH2D * m_h_W_E_vs_lhood;

      TH2D * m_h_b_pT_vs_lhood;
      TH2D * m_h_b_eta_vs_lhood;
      TH2D * m_h_b_phi_vs_lhood;
      TH2D * m_h_b_E_vs_lhood;

      TH2D * m_h_q1_pT_vs_lhood;
      TH2D * m_h_q1_eta_vs_lhood;
      TH2D * m_h_q1_phi_vs_lhood;
      TH2D * m_h_q1_E_vs_lhood;

      TH2D * m_h_q2_pT_vs_lhood;
      TH2D * m_h_q2_eta_vs_lhood;
      TH2D * m_h_q2_phi_vs_lhood;
      TH2D * m_h_q2_E_vs_lhood;

      TH2D * m_h_theta_bW_vs_lhood;
      TH2D * m_h_theta_qq_vs_lhood;

      TH2D * m_h_Rqq_vs_lhood;
      TH2D * m_h_RbW_vs_lhood;

      TH2D * m_h_t_truth_dR_vs_lhood;
      TH2D * m_h_t_truth_dPt_vs_lhood;
      TH2D * m_h_t_truth_dEta_vs_lhood;
      TH2D * m_h_t_truth_dPhi_vs_lhood;

      TH2D * m_h_b_truth_dR_vs_lhood;
      TH2D * m_h_b_truth_dPt_vs_lhood;
      TH2D * m_h_b_truth_dEta_vs_lhood;
      TH2D * m_h_b_truth_dPhi_vs_lhood;

      TH2D * m_h_W_truth_dR_vs_lhood;
      TH2D * m_h_W_truth_dPt_vs_lhood;
      TH2D * m_h_W_truth_dEta_vs_lhood;
      TH2D * m_h_W_truth_dPhi_vs_lhood;

      TLorentzVector _fjet;
      const TLorentzVector * _t_truth;
      const TLorentzVector * _b_truth;
      const TLorentzVector * _W_truth;
      const TLorentzVector * _q1_truth;
      const TLorentzVector * _q2_truth;

      TLorentzVector _b;
      TLorentzVector _q1;
      TLorentzVector _q2;
      TLorentzVector _W;
      TLorentzVector _t;

      const fastjet::ClusterSequence * _cs;
      fastjet::PseudoJet _jet;
      
};
// ---------------------------------------------------------

#endif

