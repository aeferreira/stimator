//
// Copyright (c) 2003 Joao Abecasis
//

#ifndef _agedo_ODESOLVER_H
#define _agedo_ODESOLVER_H

#include "ODESystem.h"

class ODESolver
{
protected:
	ODESystem& ODEs;
public:
	ODESolver(ODESystem& sys) : ODEs(sys)
	{
	}

	virtual int Solve(const double &time) = 0;
};

#endif
