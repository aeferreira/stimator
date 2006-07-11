
//
// Copyright (c) 2006 Antonio Ferreira
//

#include <ctime>
#include <string>
#include <ga/ga.h>
#include "DE/GADifferentialEvolution.h"
#include "ODESystem/ODESolver.h"
#include "ParametersLimits.h"
#include "Resultset.h"
#include "Operadores.h"

#include "LSODA/LivermoreSolver.h"

    int verbose_level =
    1
    ;
    int nTCs =
    1
    ;
    unsigned npars =
    4
    ;
    int nvars =
    2
    ;
    long individuals =
    80
    ;
    long generations =
    200
    ;
    char * filelist [] = {
"C:\\Documents and Settings\\tonho\\My Documents\\stimatordev\\test_at\\data1.txt"
    };

    char * filenameslist [] = {
"data1.txt"
    };

    char * parameternames [] = {
           "vmax2","vmax1","km1","km2"
    };

    char * varnames [] = {
           "B","A"
    };

    double constraints_min [] = {
           0.001,0.001,0.1,0.1
    };

    double constraints_max [] = {
           0.1,0.4,2.0,1.0
    };


class stimatorODEsystem : public ODESystem
{

	public:

	static double scale; // time scale

	stimatorODEsystem(): ODESystem(npars, nvars){}

	protected:
	virtual void
  CalculateDerivatives(const double & solution_time,
                       const double *variables, double *derivatives) const
	{

double p = 0;
if(solution_time*scale >= 10) p = 1;
double v_r1 = m_Parameters[1] * variables[1] / (m_Parameters[2] + variables[1]);
double v_r2 = p * m_Parameters[0] * variables[0] / (m_Parameters[3] + variables[0]);
derivatives[0] = +1 * v_r1 -1 * v_r2 ;
derivatives[1] = -1 * v_r1 ;


      for (int i = 0; i < nvars; ++i) derivatives[i] *= scale;

	}
};

double stimatorODEsystem::scale = 1.0;

GABoolean
GATerminateUponGenOrVarRange (GAGeneticAlgorithm & ga){
  GABoolean val = gaTrue;
  if (ga.generation() >= ga.nGenerations())
    return val;
  double ming  = ga.statistics().current(GAStatistics::Minimum);
  double maxg  = ga.statistics().current(GAStatistics::Maximum);
  double meang = ga.statistics().current(GAStatistics::Mean);

  val = gaFalse;

  if(maxg == 0) return val;
  double ratio = ((maxg-ming)/meang);
  if(ratio <= 0.000001)
      val = gaTrue;
  else
      val = gaFalse;
  return val;
}


int main(){

	//double input;

  cout << "AGEDO 1.1 Copyright (c) 2005 J.Abecasis and A.Ferreira\n\n";

  Limits lim(npars, nTCs);

  for (std::size_t i = 0; i < npars; ++i){
        lim.Min[i] = constraints_min [i];
        lim.Max[i] = constraints_max [i];
        lim.AbsoluteMin[i] = lim.AbsoluteMax[i] = true;
        if (verbose_level>=4){
            cout << "Constraints for parameter " << parameternames [i]
            << ": \n\tmin : " << lim.Min[i]
            << "\n\tmax : "   << lim.Max[i] << endl;
            cout.flush();
        }
  }

  //
  // Read result sets from file
  //
  unsigned nlines;
  for (int i = 0; i < nTCs; ++i){
      lim.ResultSets[i].nVariables = nvars + 1;
      cout << "reading time course " << filenameslist[i] << "..." ;
      //string filename(argv[i+1]);
      nlines = lim.ResultSets[i].ReadFromFile(filelist[i]);
      //nlines = lim.ResultSets[i].ReadFromFile(filename.c_str());
      cout << "\n     OK, " << nlines +1 << " lines and "
      << lim.ResultSets[i].Results.size() << " data read." << endl;
      lim.ResultSets[i].computeStatistics();
  }

	// compute timescale
	double scale = 0.0;
	for (unsigned i = 0; i < lim.ResultSets.size(); ++i){
		//Get last time of ResultSet
		double value = lim.ResultSets[i].Results[lim.ResultSets[i].Results.size() -
                   lim.ResultSets[i].nVariables].Value;
		if (value > scale) scale = value;
	}

	cout << "\nLargest time = " << scale << endl;
	stimatorODEsystem::scale = 1.0/(10 *scale);

	GA1DArrayGenome<double> genome(lim.nParameters,
          AGEDO::Objective<stimatorODEsystem, LivermoreSolver>, &lim);

	genome.initializer(AGEDO::Initialize);
	// function step was rewritten: these are not necessary in Differential evolution
  //genome.mutator(AGEDO::Mutate);
	//genome.crossover(AGEDO::SexualCrossover);

	GADifferentialEvolution ga(genome);
	ga.populationSize(individuals);
	ga.nGenerations(generations);

  cout << "\nDifferential evolution with " << ga.populationSize() <<
  " individuals and " << ga.nGenerations() << " generations \n"; cout.flush();

	ga.minimize();
	ga.flushFrequency(10);
	ga.scoreFilename("");
	ga.terminator(GATerminateUponGenOrVarRange);
	ga.selectScores(GAStatistics::AllScores);

  //
  // Generate initial solution
  //
	cout << "\nInitializing...\n";cout.flush();
  ga.initialize(std::time(NULL));
  cout << "OK\n\n";cout.flush();

    //
    // Evolve
    //

    cout << "Evolving... " << endl; cout.flush();

    while (!ga.done()){
        ga.step();
        double range  = ga.statistics().current(GAStatistics::Maximum)-
                        ga.statistics().current(GAStatistics::Minimum);
        //double meang = ga.statistics().current(GAStatistics::Mean);

        if (verbose_level >= 1){
          cout << "generation " << ga.statistics().generation()
          << ", best score = " << ga.statistics().minEver() << ",   range = " << range << "\n";
          cout.flush();
        }
        else{
            cout << "generation " << ga.statistics().generation()
            << "\n";cout.flush();
        }

        if (verbose_level >= 3){
          GA1DArrayGenome<double> const & best_individual =
              (GA1DArrayGenome<double> const &) ga.statistics().bestIndividual();

          for (size_t i = 0; i < lim.nParameters; ++i)
          {
              double value = best_individual.gene(i);
              cout << parameternames[i] << "=" << value << "\t";
          }

          cout << "\n";cout.flush();
        }
    }

  if (verbose_level >= 4) cout << ga.statistics() << endl;
	cout << "Best result:  ";
  GA1DArrayGenome<double> const & best_individual =
              (GA1DArrayGenome<double> const &) ga.statistics().bestIndividual();
  for (size_t i = 0; i < npars; ++i){
       double value = best_individual.gene(i);
       cout << parameternames[i] << " = " << value << "\t";
  }
	cout << endl;
  cout << "DONE!";

  //
  // Output best solution and associated time courses to 'best.dat'
  //
  cout << "\nGenerating best.dat...";cout.flush();
  ofstream best_file("best.dat");

  best_file.precision(9);

  best_file << "Score :\t" << best_individual.score() << endl;
  best_file << "Parameters:" << endl;
  best_file << "Name\tValue\tSE\n";

  for (size_t i = 0; i < lim.nParameters; ++i){
          double value = best_individual.gene(i);
          best_file << parameternames[i] << "\t" << value << "\tN/A\n";
  }

  best_file << endl;
  best_file << "Time courses used:\n";

  for (std::size_t i = 0; i < lim.ResultSets.size(); ++i)
      best_file << filenameslist [i] << '\n';

  best_file << endl;

  for (std::size_t i = 0; i < lim.ResultSets.size(); ++i){
      stimatorODEsystem model = stimatorODEsystem();

      // set parameters for solver
      for (unsigned j = 0; j < npars; ++j)
          model.SetParameter(j, best_individual.gene(j));

      // header
      best_file << "Time course " << filenameslist [i] << '\n';
      best_file << "time";
      for (int j = 0; j < nvars; ++j){
          best_file << '\t' << varnames[j] << "\tpred_" << varnames[j];
      }
      best_file << '\n';

      model.SetTime(lim.ResultSets[i].Results[0].Value / stimatorODEsystem::scale);
      best_file <<  lim.ResultSets[i].Results[0].Value;

      // initial values
      best_file.setf(std::ios_base::scientific);

      for (int j = 0; j < nvars; ++j){
          model.SetVariable(j, lim.ResultSets[i].Results[j + 1].Value);
          best_file << '\t' << lim.ResultSets[i].Results[j + 1].Value
          << '\t' << model.GetVariable(j);
      }
      best_file << '\n';

      LivermoreSolver solver(model);

      for (std::size_t j = nvars + 1;
           j < lim.ResultSets[i].Results.size();
           j += nvars + 1
      ){
          solver.Solve(lim.ResultSets[i].Results[j].Value/stimatorODEsystem::scale);
          best_file.unsetf(std::ios_base::scientific);
          best_file << lim.ResultSets[i].Results[j].Value;

          best_file.setf(std::ios_base::scientific);
          for (int k = 1; k < nvars + 1; ++k){
              if (!lim.ResultSets[i].Results[j + k].NotANumber)
                 best_file << '\t' << lim.ResultSets[i].Results[j + k].Value;
              else best_file << "\tNAN";
              best_file << '\t' << model.GetVariable(k-1);
          }
          best_file << '\n';
      }
      best_file << '\n';
  }
  best_file.close();
  return 0;
}

