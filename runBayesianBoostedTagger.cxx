// ***************************************************************
// This file was created using the CreateProject.sh script
// for project BayesianBoostedTagger.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

#include "BayesianBoostedTagger.h"

int main()
{

   // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   // open log file
   BCLog::OpenLog("log.txt");
   BCLog::SetLogLevel(BCLog::detail);

   // create new BayesianBoostedTagger object
   BayesianBoostedTagger * m = new BayesianBoostedTagger();

   BCLog::OutSummary("Test model created");

   // create a new summary tool object
   BCSummaryTool * summary = new BCSummaryTool(m);

   // perform your analysis here

   // normalize the posterior, i.e. integrate posterior
   // over the full parameter space
//   m->Normalize();

   // run MCMC and marginalize posterior wrt. all parameters
   // and all combinations of two parameters
//   m->MarginalizeAll();

   // run mode finding; by default using Minuit
//   m->FindMode();

   // if MCMC was run before (MarginalizeAll()) it is
   // possible to use the mode found by MCMC as
   // starting point of Minuit minimization
//   m->FindMode( m->GetBestFitParameters() );

   // draw all marginalized distributions into a PostScript file
//   m->PrintAllMarginalized("BayesianBoostedTagger_plots.ps");

   // print all summary plots
//   summary->PrintParameterPlot("BayesianBoostedTagger_parameters.eps");
//   summary->PrintCorrelationPlot("BayesianBoostedTagger_correlation.eps");
//   summary->PrintKnowledgeUpdatePlots("BayesianBoostedTagger_update.ps");

   // calculate p-value
//   m->CalculatePValue( m->GetBestFitParameters() );

   // print results of the analysis into a text file
//   m->PrintResults("BayesianBoostedTagger_results.txt");

   delete m;
   delete summary;

   BCLog::OutSummary("Test program ran successfully");
   BCLog::OutSummary("Exiting");

   // close log file
   BCLog::CloseLog();

   return 0;

}

