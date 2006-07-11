//
// Copyright (c) 2003 Joao Abecasis
//

#ifndef _agedo_GA_DIFFERENTIAL_EVOLUTION_H
#define _agedo_GA_DIFFERENTIAL_EVOLUTION_H

#include <ga/GABaseGA.h>
#include <ga/garandom.h>

class GADifferentialEvolution : public GAGeneticAlgorithm
{
public:
	GADefineIdentity("GADifferentialEvolution", 200);

	typedef int (*DifferentialCrossover)(GAGenome &p1, GAGenome &p2, GAGenome &p3, GAGenome &p4, GAGenome *c1);
	typedef double (*KGetter)(GADifferentialEvolution*);

	static int Rand1Bin(GAGenome &, GAGenome &, GAGenome &, GAGenome &, GAGenome *);
	static int CurrentToRand1(GAGenome &, GAGenome &, GAGenome &, GAGenome &, GAGenome *);
	static int Rand1Exp(GAGenome &, GAGenome &, GAGenome &, GAGenome &, GAGenome *);
	static int CurrentToRand1Bin(GAGenome &, GAGenome &, GAGenome &, GAGenome &, GAGenome *);

	static double StaticK(GADifferentialEvolution* ga) { return ga->K; }
	static double RandomK(GADifferentialEvolution*) { return GARandomDouble(); }

	GADifferentialEvolution(const GAGenome &g, DifferentialCrossover dmut = 0, double _F = 0.5, double _K = 0.5, KGetter _GetK = 0) : GAGeneticAlgorithm(g), pDiffCrossover(dmut), GetK(_GetK), K(_K), F(_F)
	{
		if (dmut == 0) pDiffCrossover = Rand1Bin;
		if (GetK == 0) GetK = StaticK;
		GetK = StaticK;
	}

	GADifferentialEvolution(const GAPopulation &p, DifferentialCrossover dmut = 0, double _F = 0.5, double _K = 0.5, KGetter _GetK = 0) : GAGeneticAlgorithm(p), pDiffCrossover(dmut), GetK(_GetK), K(_K), F(_F)
	{
		if (dmut == 0) pDiffCrossover = Rand1Bin;
		if (GetK == 0) GetK = StaticK;
		GetK = StaticK;
	}

	GADifferentialEvolution(const GADifferentialEvolution &ga) : GAGeneticAlgorithm(ga){ copy(ga); }

	virtual void copy(const GAGeneticAlgorithm &ga) {
		GAGeneticAlgorithm::copy(ga);
		const GADifferentialEvolution &gade = DYN_CAST(const GADifferentialEvolution&, ga);
		pDiffCrossover = gade.pDiffCrossover;
		GetK = gade.GetK;
		K = gade.K;
		F = gade.F;
	}

	GADifferentialEvolution& operator=(const GADifferentialEvolution &ga) {
		if (&ga != this) copy(ga);
		return *this;
	}

	virtual ~GADifferentialEvolution() { }

	virtual void initialize(unsigned int seed=0) { 
		GARandomSeed(seed);
		pop->initialize();
		pop->evaluate(gaTrue);
		stats.reset(*pop);
	}

	virtual void step(void);
	GADifferentialEvolution& operator++() { step(); return *this; }

//	static GAParameterList& registerDefaultParameters(GAParameterList&);
//	virtual int setptr(const char* name, const void* value);
//	virtual int get(const char* name, void* value) const;

public:
	DifferentialCrossover crossover() const { return pDiffCrossover; }
	DifferentialCrossover crossover(DifferentialCrossover pCross) { return pDiffCrossover = pCross; }

	KGetter GetK;	
	double K;
	double F;

protected:
	DifferentialCrossover pDiffCrossover;

};

#endif
