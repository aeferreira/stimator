//
// Copyright (c) 2003 Joao Abecasis
//

#ifndef _agedo_NUCLEAR_DECAY_H
#define _agedo_NUCLEAR_DECAY_H

#include "ODESystem.h"

class NuclearDecay : public ODESystem
{
	// Parâmetros
	//  1.	k1
	//  2.	k2

	// Variáveis dinâmicas
	//  1.	x
	//  2.	y
	//  3.	z

	public:
	// constructor
	NuclearDecay(double _k1 = 0., double _k2 = 0.) : ODESystem(2, 3)
	{
		m_Parameters[0] = _k1;
		m_Parameters[1] = _k2;
	}

	protected:
	// ODEs
	virtual void CalculateDerivatives(const double &/*time*/, const double *vars, double *derivs) const
	// Two step nuclear decay:
	//   x -> y (k1)
	//   y -> z (k2)
	{
		derivs[0] = - m_Parameters[0] * vars[0];
		derivs[2] = m_Parameters[1] * vars[1];
		derivs[1] = - derivs[0] - derivs[2];
	}
};

#endif
