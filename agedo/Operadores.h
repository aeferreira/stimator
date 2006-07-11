//
// Copyright (c) 2004 Joao Abecasis / António Ferreira
//

#include <cmath>
#include <limits>

#include "ODESystem/ODESolver.h"
#include "ParametersLimits.h"

namespace AGEDO {

	using namespace std;

	void Initialize(GAGenome& g)
	{
		GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)g;
		Limits * lim = (Limits*) g.userData();

		for (unsigned int i = 0; i < (unsigned int)genome.length(); ++i)
		{
			genome.gene(i, GARandomDouble(lim->Min[i], lim->Max[i]));
		}
		return;
	}

	int Mutate(GAGenome &g, float pmut)
	{
		if (pmut <= 0) return 0;

		int nmut = 0;

		GA1DArrayGenome<double> &child((GA1DArrayGenome<double>&) g);
		Limits * lim = (Limits*) g.userData();

		for (unsigned int i = 0; i < (unsigned int)child.length(); i++)
		{
			if (GAFlipCoin(pmut))
			{
				if (GAFlipCoin(0.5))
				{
					child.gene(i, GARandomDouble(child.gene(i), lim->Max[i]));
				}
				else
				{
					child.gene(i, GARandomDouble(lim->Min[i], child.gene(i)));
				}
				nmut++;
			}
		}

		return nmut;
	}

	int SexualCrossover(const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2)
	{
		int nCross = 0;

		const GA1DArrayGenome<double> &mom((const GA1DArrayGenome<double>&) p1);
		const GA1DArrayGenome<double> &dad((const GA1DArrayGenome<double>&) p2);
		if(c1)
		{
			GA1DArrayGenome<double> &sis((GA1DArrayGenome<double>&) *c1);
			for (unsigned int i = 0; i < (unsigned int)mom.length(); i++)
			{
				if (GAFlipCoin(0.7))
				{
					sis.gene(i, mom.gene(i) + GARandomDouble() * (dad.gene(i) - mom.gene(i)));
					nCross++;
				}
				else
				{
					sis.gene(i, mom.gene(i));
				}
			}
		}

		if(c2)
		{
			GA1DArrayGenome<double> &bro((GA1DArrayGenome<double>&) *c2);
			for (unsigned int i = 0; i < (unsigned int)mom.length(); i++)
			{
				if (GAFlipCoin(0.7))
				{
					bro.gene(i,mom.gene(i) + GARandomDouble() * (dad.gene(i) - mom.gene(i)));
					nCross++;
				}
				else
				{
					bro.gene(i, dad.gene(i));
				}
			}
		}
		return nCross;
	}

	template<class DerivedFromODESystem, class DerivedFromODESolver>
	float Objective(GAGenome& g)
	{
		unsigned int i, j, k;
		GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)g;
		Limits * lim = (Limits*) g.userData();

		DerivedFromODESystem odes;

		for (i = 0; i < (unsigned int)genome.length(); ++i)
		{
			odes.SetParameter(i, genome.gene(i));
		}

		double * runResult = new double[lim->ResultSets.size()];
		double sum, square;

		for (i = 0; i < lim->ResultSets.size(); ++i)
		{
			sum = 0.;
			if (lim->ResultSets[i].discard) {
				runResult[i] = 0.0;
				continue;
			}
			// Set initial values on odes
			odes.SetTime(lim->ResultSets[i].Results[0].Value / DerivedFromODESystem::scale);
			for (j = 0; j < odes.nVariables(); ++j)
			{
				// missing i.v. are set to zero
        if (!lim->ResultSets[i].Results[j + 1].NotANumber)
            odes.SetVariable(j, lim->ResultSets[i].Results[j + 1].Value);
        else
            odes.SetVariable(j, 0.0);
			}

			{// Better to get a new solver every time...
				DerivedFromODESolver solver(odes);
				for (j = odes.nVariables() + 1; j < lim->ResultSets[i].Results.size(); j += odes.nVariables() + 1)
				{
					//Iterate through remaining values
					if (solver.Solve(lim->ResultSets[i].Results[j].Value / DerivedFromODESystem::scale)  < 0)
					{
						// TODO: Remove predefined penalty
            // runResult[i] = sum + 5000.;
                        sum = 5000.;
						break;
					}
					for (k = 0; k < odes.nVariables(); ++k)
					{
						if (!lim->ResultSets[i].Results[j + k + 1].NotANumber)
						{
							double expvalue = lim->ResultSets[i].Results[j + k + 1].Value;
							square = odes.GetVariable(k) - expvalue;
							square *= square;
							sum += square;
						}
					}
				}
			}
			runResult[i] = sum;
		}

		sum = 0.0;
		for (i = 0; i < lim->ResultSets.size(); ++i)
			sum += runResult[i];

		float ret_val;
		ret_val = min(sum, (double)numeric_limits<float>::max());
		ret_val = max(sum, (double)numeric_limits<float>::min());
    delete[] runResult;
		return ret_val;
	}


    template<class DerivedFromODESystem, class DerivedFromODESolver>
    float ObjectiveM(GAGenome& g)
    {
        using std::min;
        using std::max;

        unsigned int i, j, k;
        GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)g;
        Limits * lim = (Limits*) g.userData();

        DerivedFromODESystem odes;

        for (i = 0; i < (unsigned int)genome.length(); ++i)
        {
            odes.SetParameter(i, genome.gene(i));
        }

        double * runResult = new double[lim->ResultSets.size()];
        double sum, square;

        for (i = 0; i < lim->ResultSets.size(); ++i)
        {
            sum = 0.;
            if (lim->ResultSets[i].discard) {
                runResult[i] = 0.0;
                continue;
            }
            // Set initial values on odes
            odes.SetTime(lim->ResultSets[i].Results[0].Value
                         / DerivedFromODESystem::scale);

            for (j = 0; j < odes.nVariables(); ++j){
        				// missing i.v. are set to zero
                if (!lim->ResultSets[i].Results[j + 1].NotANumber)
                    odes.SetVariable(j, lim->ResultSets[i].Results[j + 1].Value);
                else
                    odes.SetVariable(j, 0.0);
            }

            {// Better to get a new solver every time...
                DerivedFromODESolver solver(odes);
                for (j = odes.nVariables() + 1;
                     j < lim->ResultSets[i].Results.size();
                     j += odes.nVariables() + 1){
                    //Iterate through remaining values
                    if (solver.Solve(lim->ResultSets[i].Results[j].Value
                                     / DerivedFromODESystem::scale) < 0
                    ){
                        // TODO: Remove predefined penalty
                        // runResult[i] = sum + 5000.;
                        sum = 5000.;
                        break;
                    }
                    for (k = 0; k < odes.nVariables(); ++k){
                        if (!lim->ResultSets[i].Results[j + k + 1].NotANumber){
                            double expvalue =
                                lim->ResultSets[i].Results[j + k + 1].Value;
                            square = odes.GetVariable(k) - expvalue;
                            square /= lim->ResultSets[i].average[k + 1];
                            square *= square;
                            sum += square;
                        }
                    }
                }
            }
            runResult[i] = sum;
        }

        sum = 0.0;
        for (i = 0; i < lim->ResultSets.size(); ++i)
            sum += runResult[i];

        float ret_val;
        ret_val = min(sum, (double)std::numeric_limits<float>::max());
        ret_val = max(sum, (double)std::numeric_limits<float>::min());
        delete[] runResult;
        return ret_val;
    }



	template<class DerivedFromODESystem, class DerivedFromODESolver>
	float ObjectiveSR(GAGenome& g)
	{
		unsigned int i, j, k;
		GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)g;
		Limits * lim = (Limits*) g.userData();

		DerivedFromODESystem odes;

		for (i = 0; i < (unsigned int)genome.length(); ++i)
		{
			odes.SetParameter(i, genome.gene(i));
		}

		double * runResult = new double[lim->ResultSets.size()];
		double sum, square;

		for (i = 0; i < lim->ResultSets.size(); ++i)
		{
			sum = 0.;
			if (lim->ResultSets[i].discard) {
				runResult[i] = 0.0;
				continue;
			}
			// Set initial values on odes
			odes.SetTime(lim->ResultSets[i].Results[0].Value);
      for (j = 0; j < odes.nVariables(); ++j){
  				// missing i.v. are set to zero
          if (!lim->ResultSets[i].Results[j + 1].NotANumber)
              odes.SetVariable(j, lim->ResultSets[i].Results[j + 1].Value);
          else
              odes.SetVariable(j, 0.0);
      }

			{// Better to get a new solver every time...
				DerivedFromODESolver solver(odes);
				for (j = odes.nVariables() + 1; j < lim->ResultSets[i].Results.size(); j += odes.nVariables() + 1)
				{
					//Iterate through remaining values
					if (solver.Solve(lim->ResultSets[i].Results[j].Value) < 0)
					{
						// TODO: Remove predefined penalty
//						runResult[i] = sum + 5000.;
                        sum = 5000.;
						break;
					}
					for (k = 0; k < odes.nVariables(); ++k)
					{
						if (!lim->ResultSets[i].Results[j + k + 1].NotANumber)
						{
							double expvalue = lim->ResultSets[i].Results[j + k + 1].Value;
							square = odes.GetVariable(k) - expvalue;
							square *= square;
							sum += square;
						}
					}
				}
			}
			runResult[i] = sum / (lim->ResultSets[i].range);
		}

		sum = 0.0;
		for (i = 0; i < lim->ResultSets.size(); ++i)
			sum += runResult[i];

		float ret_val;
		ret_val = min(sum, (double)numeric_limits<float>::max());
		ret_val = max(sum, (double)numeric_limits<float>::min());
        delete[] runResult;
		return ret_val;
	}

}

