//
// Copyright (c) 2003 Joao Abecasis
//

#include "GADifferentialEvolution.h"
#include "../ParametersLimits.h"
#include <ga/GA1DArrayGenome.h>

#if 0
GADifferentialEvolution::GADifferentialEvolution(const GAGenome& g) : GAGeneticAlgorithm(g)
{
	// Set Default Parameters
	// Add Parameters to list
	//  params.add(FullName, ShortName, Type, DefaultValue);
}

GADifferentialEvolution::GADifferentialEvolution(const GAPopulation& p) : GAGeneticAlgorithm(p)
{
	// Set Default Parameters
	// Add Parameters to list
	//  params.add(FullName, ShortName, Type, DefaultValue);
}

GADifferentialEvolution::GADifferentialEvolution(const GADifferentialEvolution& ga) : GAGeneticAlgorithm(ga)
{
	copy(ga);
}
#endif

void GADifferentialEvolution::step(void)
{
	int i, r1, r2, r3, k;

	// Clone an existing genome for temp work...
	GAGenome *g = pop->individual(0).clone();

	for (i = 0; i < pop->size(); ++i)
	{
		do { r1 = (GARandomInt(0, pop->size() - 1)); } while ( r1 == i );
		do { r2 = (GARandomInt(0, pop->size() - 1)); } while ( r2 == i || r2 == r1 );
		do { r3 = (GARandomInt(0, pop->size() - 1)); } while ( r3 == i || r3 == r1 || r3 == r2 );

		// Crossover/Mutate
		if ((k = (*pDiffCrossover)(pop->individual(i), pop->individual(r1), pop->individual(r2), pop->individual(r3), g)) != 0)
		{
			stats.nummut += k;
//           cout << "g = ";             
//           for (int ii = 0; ii < ((GA1DArrayGenome<double>*)g)->length() ; ++ii)
//                   cout << ((GA1DArrayGenome<double>*)g)->gene(ii) << " ";
//           cout << endl;
//           cout << i << " = ";             
//           for (int ii = 0; ii < ((GA1DArrayGenome<double>*)g)->length() ; ++ii)
//                   cout << ((GA1DArrayGenome<double>*)g)->gene(ii) << " ";
//           cout << endl;
			if (g->evaluate(gaTrue) != pop->individual(i).score())
			{
				if (((minimaxi() == GAGeneticAlgorithm::MINIMIZE) && (g->score() < pop->individual(i).score())) || ((minimaxi() == GAGeneticAlgorithm::MAXIMIZE) && (g->score() > pop->individual(i).score())))
				{
					g = pop->replace(g, i);
					++stats.numrep;
				}
			}
		}
		++stats.numcro;
	}
	stats.numeval += pop->size();

	stats.update(*pop);

	delete g;
	
	return;
}

int GADifferentialEvolution::Rand1Bin(GAGenome &p1, GAGenome &p2, GAGenome &p3, GAGenome &p4, GAGenome *c1)
{// This is the Original version
	if (c1 == 0) return 0;

	GADifferentialEvolution* ga = (GADifferentialEvolution*) p1.geneticAlgorithm();
	if (ga == 0) return 0;

	GA1DArrayGenome<double>& mom = (GA1DArrayGenome<double>&) p1;
	GA1DArrayGenome<double>& dad1 = (GA1DArrayGenome<double>&) p2;
	GA1DArrayGenome<double>& dad2 = (GA1DArrayGenome<double>&) p3;
	GA1DArrayGenome<double>& dad3 = (GA1DArrayGenome<double>&) p4;
	GA1DArrayGenome<double>& bro = *((GA1DArrayGenome<double>*)c1);

	Limits * lim = (Limits*) mom.userData();
	double temp;

	int j = GARandomInt(0, mom.length() - 1), l = 1;
	for (int k = mom.length(); k > 0; --k)
	{
		if ((GAFlipCoin(ga->pCrossover())))
		{
			temp = dad3.gene(j) + ga->F * (dad1.gene(j) - dad2.gene(j));
			if (lim->AbsoluteMax[j]) temp = min(temp, lim->Max[j]);
			if (lim->AbsoluteMin[j]) temp = max(temp, lim->Min[j]);
			bro.gene(j, temp);
			++l;
		}
		else
		{
			bro.gene(j, mom.gene(j));
		}
		++j %= mom.length();
	}
	temp = dad3.gene(j) + ga->F * (dad1.gene(j) - dad2.gene(j));
	if (lim->AbsoluteMax[j]) temp = min(temp, lim->Max[j]);
	if (lim->AbsoluteMin[j]) temp = max(temp, lim->Min[j]);
	bro.gene(j, temp);

	return l;
}

int GADifferentialEvolution::CurrentToRand1(GAGenome &p1, GAGenome &p2, GAGenome &p3, GAGenome &p4, GAGenome *c1)
{// Every gene Mutates!
	if (c1 == 0) return 0;

	GADifferentialEvolution* ga = (GADifferentialEvolution*) p1.geneticAlgorithm();
	if (ga == 0) return 0;

	GA1DArrayGenome<double>& mom = (GA1DArrayGenome<double>&) p1;
	GA1DArrayGenome<double>& dad1 = (GA1DArrayGenome<double>&) p2;
	GA1DArrayGenome<double>& dad2 = (GA1DArrayGenome<double>&) p3;
	GA1DArrayGenome<double>& dad3 = (GA1DArrayGenome<double>&) p4;
	GA1DArrayGenome<double>& bro = *((GA1DArrayGenome<double>*)c1);

	Limits * lim = (Limits*) mom.userData();
	double temp;

	for (int j = mom.length() - 1; j >= 0; --j)
	{
		temp = mom.gene(j) + ga->GetK(ga) * (dad3.gene(j) - mom.gene(j)) + ga->F * (dad1.gene(j) - dad2.gene(j));
		if (lim->AbsoluteMax[j]) temp = min(temp, lim->Max[j]);
		if (lim->AbsoluteMin[j]) temp = max(temp, lim->Min[j]);
		bro.gene(j, temp);
	}

	return mom.length();
}

int GADifferentialEvolution::Rand1Exp(GAGenome &p1, GAGenome &p2, GAGenome &p3, GAGenome &p4, GAGenome *c1)
{// Once crossover begins, there's no stopping it
	if (c1 == 0) return 0;

	GADifferentialEvolution* ga = (GADifferentialEvolution*) p1.geneticAlgorithm();
	if (ga == 0) return 0;

	GA1DArrayGenome<double>& mom = (GA1DArrayGenome<double>&) p1;
	GA1DArrayGenome<double>& dad1 = (GA1DArrayGenome<double>&) p2;
	GA1DArrayGenome<double>& dad2 = (GA1DArrayGenome<double>&) p3;
	GA1DArrayGenome<double>& dad3 = (GA1DArrayGenome<double>&) p4;
	GA1DArrayGenome<double>& bro = *((GA1DArrayGenome<double>*)c1);

	Limits * lim = (Limits*) mom.userData();
	double temp;

	int j = GARandomInt(0, mom.length() - 1), l = 1;
	bool flag = false;
	for (int k = mom.length(); k > 0; --k)
	{
		if (GAFlipCoin(ga->pCrossover())) flag = true;
		if (flag)
		{
			temp = dad3.gene(j) + ga->F * (dad1.gene(j) - dad2.gene(j));
			if (lim->AbsoluteMax[j]) temp = min(temp, lim->Max[j]);
			if (lim->AbsoluteMin[j]) temp = max(temp, lim->Min[j]);
			bro.gene(j, temp);
			++l;
		}
		else
			bro.gene(j, mom.gene(j));
		++j %= mom.length();
	}

	temp = dad3.gene(j) + ga->F * (dad1.gene(j) - dad2.gene(j));
	if (lim->AbsoluteMax[j]) temp = min(temp, lim->Max[j]);
	if (lim->AbsoluteMin[j]) temp = max(temp, lim->Min[j]);
	bro.gene(j, temp);

	return l;
}

int GADifferentialEvolution::CurrentToRand1Bin(GAGenome &p1, GAGenome &p2, GAGenome &p3, GAGenome &p4, GAGenome *c1)
{
	if (c1 == 0) return 0;

	GADifferentialEvolution* ga = (GADifferentialEvolution*) p1.geneticAlgorithm();
	if (ga == 0) return 0;

	GA1DArrayGenome<double>& mom = (GA1DArrayGenome<double>&) p1;
	GA1DArrayGenome<double>& dad1 = (GA1DArrayGenome<double>&) p2;
	GA1DArrayGenome<double>& dad2 = (GA1DArrayGenome<double>&) p3;
	GA1DArrayGenome<double>& dad3 = (GA1DArrayGenome<double>&) p4;
	GA1DArrayGenome<double>& bro = *((GA1DArrayGenome<double>*)c1);

	Limits * lim = (Limits*) mom.userData();
	double temp;

	int j = GARandomInt(0, mom.length() - 1), l = 1;
	for (int k = mom.length(); k > 0; --k)
	{
		if (GAFlipCoin(ga->pCrossover()))
		{
			temp = mom.gene(j) + ga->GetK(ga) * (dad3.gene(j) - mom.gene(j)) + ga->F * (dad1.gene(j) - dad2.gene(j));
			if (lim->AbsoluteMax[j]) temp = min(temp, lim->Max[j]);
			if (lim->AbsoluteMin[j]) temp = max(temp, lim->Min[j]);
			bro.gene(j, temp);
			++l;
		}
		else
			bro.gene(j, mom.gene(j));
		++j %= mom.length();
	}

	temp = mom.gene(j) + ga->GetK(ga) * (dad3.gene(j) - mom.gene(j)) + ga->F * (dad1.gene(j) - dad2.gene(j));
	if (lim->AbsoluteMax[j]) temp = min(temp, lim->Max[j]);
	if (lim->AbsoluteMin[j]) temp = max(temp, lim->Min[j]);
	bro.gene(j, temp);

	return l;
}

