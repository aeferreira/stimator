// $Header: //194.117.40.12/c\044/Documents\040and\040Settings/CVSNT/Repositories/Bioquimica\040Teorica/AGEDO/Source/ga/GADCrowdingGA.h,v 1.1 2003/09/05 16:46:14 Administrator Exp $
/* ----------------------------------------------------------------------------
  dcrowdingga.h
  mbwall 29mar99
  Copyright (c) 1999 Matthew Wall, all rights reserved

  Header file for the steady-state genetic algorithm class.
---------------------------------------------------------------------------- */
#ifndef _ga_deterministic_crowding_ga_h_
#define _ga_deterministic_crowding_ga_h_

#include <ga/GABaseGA.h>

class GADCrowdingGA : public GAGeneticAlgorithm {
public:
  GADefineIdentity("GADeterministicCrowdingGA", 241);

  GADCrowdingGA(const GAGenome& g) : GAGeneticAlgorithm(g) {}
  virtual ~GADCrowdingGA() {}

  virtual void initialize(unsigned int seed=0);
  virtual void step();
  GADCrowdingGA& operator++() { step(); return *this; }
};

#endif
