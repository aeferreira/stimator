//
// Copyright (c) 2003 Joao Abecasis
//

#ifndef _agedo_GLIOXALASES_H
#define _agedo_GLIOXALASES_H

#include "ODESystem.h"

class Glioxalases : public ODESystem
{
	// Modelo cinético da via dos glioxalases (pág 118)
	//    	           k1
	//  1.	MG + GSH <====> HTA
	//    	           k2
	//
	//    	     GlxI
	//  2.	HTA ------> SDLGSH
	// 
	//    	        GlxII
	//  3.  SDLGSH -------> L + GSH

	// Reagentes
	//  1.	MG	Metilglioxal
	//  2.	GSH	Glutationo (L-gama-glutamilo-L-cisteinilo-glicina)
	//  3.	HTA	Hemitioacetal
	//  4.	SDLGSH	S-D-Lactoilglutationo
	//  5.	L	D-Lactato

	// Enzimas envolvidos
	//  1.	GlxI	Glioxalase I (S-D-lactoilglutationo:metilglioxal liase, EC 4.4.1.5)
	//  	Cinética HMM irreversível
	//  2.	GlxII	Glioxalase II (S-2-hidroxiacilglutationo hidrolase, EC 3.1.2.6)
	//  	Cinética HMM irreversível

	// Parâmetros
	//  1.	k1
	//  2.	k2
	//  3.	Vglx1
	//  4.	Kglx1
	//  5.	Vglx2
	//  6.	Kglx2

	// Referências
	//  1.	Ana Martins (2000), O sistema dos Glioxalases e Metabolismo do 
	//  	Metilglioxal em Saccharomyces cerevisiae, Tese de doutoramento, 
	//  	Universidade de Lisboa, Lisboa

	public:
	// constructor
	Glioxalases(double _k1 = 0., double _k2 = 0., double _Vglx1 = 0., double _Kglx1 = 0., double _Vglx2 = 0., double _Kglx2 = 0.) : ODESystem(6, 5)
	{
		m_Parameters[0] = _k1;
		m_Parameters[1] = _k2;
		m_Parameters[2] = _Vglx1;
		m_Parameters[3] = _Kglx1;
		m_Parameters[4] = _Vglx2;
		m_Parameters[5] = _Kglx2;
	}

	protected:
	// ODEs
	virtual void CalculateDerivatives(const double &/* time */, const double *vars, double *derivs) const
	{
		derivs[0] = - m_Parameters[0] * vars[0] * vars[1] + m_Parameters[1] * vars[2];
		derivs[4] = m_Parameters[4] * vars[3] / (m_Parameters[5] + vars[3]);
		derivs[1] = derivs[0] + derivs[4];
		derivs[2] = - derivs[0] - m_Parameters[2] * vars[2] / (m_Parameters[3] + vars[2]);
		derivs[3] = m_Parameters[2] * vars[2] / (m_Parameters[3] + vars[2]) - derivs[4];
	}
};

#endif
