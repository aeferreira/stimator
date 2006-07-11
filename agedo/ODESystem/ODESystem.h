//
// Copyright (c) 2003 Joao Abecasis
//

#ifndef _agedo_ODE_SYSTEM_H
#define _agedo_ODE_SYSTEM_H

#include <cassert>
#include <cmath>

class ODESystem
{
	protected:
	unsigned int m_nParameters;			// Number of Static Variables
	unsigned int m_nVariables;			// Number of Dynamic Variables
	double m_curTime;			// Current Time
	double * m_Parameters;			// Static Variables
	double * m_curVariables;		// Current Value of Dynamic Variables
	mutable double * m_curDerivatives;	// Current Value of Variables' Derivatives
	mutable bool m_bDerivativesCalculated;

	public:
	ODESystem(unsigned int nParams = 0, unsigned int nVars = 0) : m_nParameters(nParams), m_nVariables(nVars), m_curTime(0.), m_Parameters(0), m_curVariables(0), m_curDerivatives(0), m_bDerivativesCalculated(false)
	{
		m_Parameters = (m_nParameters == 0) ? 0 : new double[m_nParameters];
		m_curVariables = (m_nVariables == 0) ? 0 : new double[m_nVariables];
		m_curDerivatives = (m_nVariables == 0) ? 0 : new double[m_nVariables];
	}

	ODESystem(const ODESystem &obj) : m_nParameters(obj.m_nParameters), m_nVariables(obj.m_nVariables), m_curTime(obj.m_curTime), m_Parameters(0), m_curVariables(0), m_curDerivatives(0), m_bDerivativesCalculated(false)
	{
		if (m_nParameters != 0)
		{
			try {
				m_Parameters = new double[m_nParameters];
			} catch (...) {
				delete[] m_Parameters;
				throw;
			}

			for (unsigned int i = 0; i < m_nParameters; ++i)
			{
				m_Parameters[i] = obj.m_Parameters[i];
			}
		}

		if (m_curVariables != 0)
		{
			try {
				m_curVariables = new double[m_nVariables];
				m_curDerivatives = new double[m_nVariables];
			} catch (...) {
				delete [] m_Parameters;
				delete [] m_curVariables;
				delete [] m_curDerivatives;
				throw;
			}

			for (unsigned int i = 0; i < m_nVariables; ++i)
			{
				m_curVariables[i] = obj.m_curVariables[i];
			}
		}
	}

	~ODESystem()
	{
		delete [] m_Parameters;
		delete [] m_curVariables;
		delete [] m_curDerivatives;
	}

	protected:
	ODESystem& operator=(const ODESystem &that)
	{
		// ODESystem shouldn't be assigned,unless you ABSOLUTELY know what you're doing, of course!
		// ODESolvers keep a reference to ODESystems and might keep some more internal information 
		// related to the system they're solving.
		// It is important that only ODESystems of the same derived class are copy assigned.
		// Anyway... we provide an implementation for when nParameters and nVariables remain 
		// constant so derived classes can concentrate on copying their own structures.
		if (this != &that)
		{
			assert(m_nParameters == that.m_nParameters);
			assert(m_nVariables == that.m_nVariables);

			m_bDerivativesCalculated = false;
			m_curTime = that.m_curTime;

			for (unsigned int i = 0; i < m_nParameters; ++i)
			{
				m_Parameters[i] = that.m_Parameters[i];
			}

			for (unsigned int i = 0; i < m_nVariables; ++i)
			{
				m_curVariables[i] = that.m_curVariables[i];
			}
		}
		return *this;
	}

	public:
	virtual unsigned int nParameters(void) const
	{
		return m_nParameters;
	}

	virtual unsigned int nVariables(void) const
	{
		return m_nVariables;
	}

	virtual void SetParameter(unsigned int index, const double &value)
	{
		assert(index < m_nParameters);
		assert(m_Parameters != 0);

		if (m_Parameters[index] != value)
		{
			m_Parameters[index] = value;
			m_bDerivativesCalculated = false;
		}
		return;
	}

	virtual double GetParameter(unsigned int index) const
	{
		assert(index < m_nParameters);
		assert(m_Parameters != 0);

		return m_Parameters[index];
	}

	virtual void SetVariable(unsigned int index, const double &value)
	{
		assert(index < m_nVariables);
		assert(m_curVariables != 0);

		if (m_curVariables[index] != value)
		{
			m_curVariables[index] = value;
			m_bDerivativesCalculated = false;
		}
		return;
	}

	virtual double GetVariable(unsigned int index) const
	{
		assert(index < m_nVariables);
		assert(m_curVariables != 0);

		return m_curVariables[index];
	}

	virtual void SetTime(const double &value)
	{
		if (m_curTime != value)
		{
			m_curTime = value;
			m_bDerivativesCalculated = false;
		}
		return;
	}

	virtual double GetTime(void) const
	{
		return m_curTime;
	}

	virtual void CalculateDerivatives(const double & time, const double *vars, double *derivs) const = 0;

//	virtual void CalculateDerivatives(void) const
//	{
//		if (!m_bDerivativesCalculated)
//		{
//			this->CalculateDerivatives(m_curTime, m_curVariables, m_curDerivatives);
//			m_bDerivativesCalculated = true;
//		}
//		return;
//	}

	virtual double GetDerivative(unsigned int index) const
	{
		assert(index < m_nVariables);
		assert(m_curVariables != 0);
		assert(m_curDerivatives != 0);

		if (!m_bDerivativesCalculated)
		{
			this->CalculateDerivatives(m_curTime, m_curVariables, m_curDerivatives);
			m_bDerivativesCalculated = true;
		}

		return m_curDerivatives[index];
	}


	virtual void Jacobian(double *)
	{
		/* TODO */
		throw "Function Not Implemented";
	}

	virtual void GetErrorWeights(const long itol, const double *const rtol, const double *const atol, const double *const ycur, double *ewt)
	{
		switch (itol)
		{
		case 1:
			for (unsigned int i = 0; i < m_nVariables; ++i) {
				ewt[i] = rtol[0] * fabs(ycur[i]) + atol[0];
			}
			break;
		case 2:
			for (unsigned int i = 0; i < m_nVariables; ++i) {
				ewt[i] = rtol[0] * fabs(ycur[i]) + atol[i];
			}
			break;
		case 3:
			for (unsigned int i = 0; i < m_nVariables; ++i) {
				ewt[i] = rtol[i] * fabs(ycur[i]) + atol[0];
			}
			break;
		case 4:
			for (unsigned int i = 0; i < m_nVariables; ++i) {
				ewt[i] = rtol[i] * fabs(ycur[i]) + atol[i];
			}
			break;
		}
		return;
	}
};

#endif
