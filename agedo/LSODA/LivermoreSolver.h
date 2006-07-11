//
// Copyright (c) 2003 Joao Abecasis
//

#ifndef _agedo_LIVERMORE_SOLVER_H
#define _agedo_LIVERMORE_SOLVER_H

#include "ErrorMessenger.h"
#include "../ODESystem/ODESolver.h"

#include <limits>
#include <cmath>

using namespace std;

class LivermoreSolver : public ODESolver
{
	LivermoreSolver &operator=(LivermoreSolver &);
	LivermoreSolver(LivermoreSolver &);

#ifdef LSODA_CALCULATE_ODE_COEFFICIENTS

	static void SetAdamsCoefficients(double * const elco, double * const tesco)
	// meth == 1
	// elco[13 * 12], tesco[3 * 12] both in column-major order
	{
	//  Set ODE integrator coefficients.
	//  
	//  AUTHOR  Hindmarsh, Alan C., (LLNL)
	//
	//  called by the integrator routine to set coefficients
	//  needed there.  The coefficients for the current method, as
	//  given by the value of METH, are set for all orders and saved.
	//  The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2.
	//  (A smaller value of the maximum order is also allowed.)
	//  DCFODE is called once at the beginning of the problem,
	//  and is not called again unless and until METH is changed.
	//
	//  The ELCO array contains the basic method coefficients.
	//  The coefficients el(i), 1 .le. i .le. nq+1, for the method of
	//  order nq are stored in ELCO(i,nq).  They are given by a generating
	//  polynomial, i.e.,
	//      l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
	//  For the implicit Adams methods, l(x) is given by
	//      dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
	//  For the BDF methods, l(x) is given by
	//      l(x) = (x+1)*(x+2)* ... *(x+nq)/K,
	//  where         K = factorial(nq)*(1 + 1/2 + ... + 1/nq).
	//
	//  The TESCO array contains test constants used for the
	//  local error test and the selection of step size and/or order.
	//  At order nq, TESCO(k,nq) is used for the selection of step
	//  size at order nq - 1 if k = 1, at order nq if k = 2, and at order
	//  nq + 1 if k = 3.

		long int i, j;
		double pc[12], d2, d3, tsign, d4;

		elco[0] = 1.;
		elco[1] = 1.;
		tesco[0] = 0.;
		tesco[1] = 2.;
		tesco[3] = 1.;
		tesco[35] = 0.;
		pc[0] = 1.;
		d3 = 1.;
		for (i = 1; i < 12; ++i)
		{
			pc[i] = 1.;
			for (j = i - 1; j > 0; --j)
			{
				pc[j] = pc[j-1] + i*pc[j];
			}
			pc[0] *= i;
			d2 = pc[0];
			d4 = pc[0]/2.;
			tsign = 1.;
			for (j = 1; j <= i; ++j)
			{
				tsign = -tsign;
				d2 += tsign*pc[j]/(j+1);
				d4 += tsign*pc[j]/(j+2);
			}
			elco[13*i] = d2*d3;
			elco[13*i + 1] = 1.;
			for (j = 1; j < i; ++j)
			{
				elco[13*i + j + 1] = d3*pc[j]/(j+1);
			}
			d3 /= (i+1);
			elco[14*i + 1] = d3;
			tesco[3*i + 1] = tesco[3*i - 1] = 1./(d3*d4);
			(i < 11) ? tesco[3*i + 3] = 1./(d4*(i + 2)) : 0;
		}
		return;
	}

	static void SetBDFCoefficients(double * const elco, double * const tesco)
	// meth == 2
	{
		long int i, j;
		double pc[6];

		pc[0] = 1.;
		for (i = 0; i < 5; ++i)
		{
			pc[i + 1] = 0.;
			for (j = i + 1; j > 0; --j)
			{
				pc[j] = pc[j-1] + (i+1)*pc[j];
			}
			tesco[3*i] = 1./pc[0];
			pc[0] *= i + 1;
			for (j = i + 1; j >= 0; --j)
			{
				elco[j+13*i] = pc[j]/pc[1];
			}
			elco[1+13*i] = 1.;
			tesco[1+3*i] = (i+2)/elco[13*i];
			tesco[2+3*i] = (i+3)/elco[13*i];
		}
		return;
	}

#else // defined(LSODA_CALCULATE_ODE_COEFFICIENTS)

	static void SetAdamsCoefficients(double* const elco, double* const tesco)
	// meth == 1
	{
		elco[0] = 1.;
		elco[1] = 1.;
		elco[13] = 1./2.;
		elco[14] = 1.;
		elco[15] = 1./2.;
		elco[26] = 5./12.;
		elco[27] = 1.;
		elco[28] = 3./4.;
		elco[29] = 1./6.;
		elco[39] = 3./8.;
		elco[40] = 1.;
		elco[41] = 11./12.;
		elco[42] = 1./3.;
		elco[43] = 1./24.;
		elco[52] = 251./720.;
		elco[53] = 1.;
		elco[54] = 25./24.;
		elco[55] = 35./72.;
		elco[56] = 5./48.;
		elco[57] = 1./120.;
		elco[65] = 95./288.;
		elco[66] = 1.;
		elco[67] = 137./120.;
		elco[68] = 5./8.;
		elco[69] = 17./96.;
		elco[70] = 1./40.;
		elco[71] = 1./720.;
		elco[78] = 19087./60480.;
		elco[79] = 1.;
		elco[80] = 49./40.;
		elco[81] = 203./270.;
		elco[82] = 49./192.;
		elco[83] = 7./144.;
		elco[84] = 7./1440.;
		elco[85] = 1./5040.;
		elco[91] = 5257./17280.;
		elco[92] = 1.;
		elco[93] = 363./280.;
		elco[94] = 469./540.;
		elco[95] = 967./2880.;
		elco[96] = 7./90.;
		elco[97] = 23./2160.;
		elco[98] = 1./1260.;
		elco[99] = 1./40320.;
		elco[104] = 1070017./3628800.;
		elco[105] = 1.;
		elco[106] = 761./560.;
		elco[107] = 29531./30240.;
		elco[108] = 267./640.;
		elco[109] = 1069./9600.;
		elco[110] = 3./160.;
		elco[111] = 13./6720.;
		elco[112] = 1./8960.;
		elco[113] = 1./362880.;
		elco[117] = 25713./89600.;
		elco[118] = 1.;
		elco[119] = 7129./5040.;
		elco[120] = 6515./6048.;
		elco[121] = 4523./9072.;
		elco[122] = 19./128.;
		elco[123] = 3013./103680.;
		elco[124] = 5./1344.;
		elco[125] = 29./96768.;
		elco[126] = 1./72576.;
		elco[127] = 1./3628800.;
		elco[130] = 26842253./95800320.;
		elco[131] = 1.;
		elco[132] = 7381./5040.;
		elco[133] = 177133./151200.;
		elco[134] = 84095./145152.;
		elco[135] = 341693./1814400.;
		elco[136] = 8591./207360.;
		elco[137] = 7513./1209600.;
		elco[138] = 121./193536.;
		elco[139] = 11./272160.;
		elco[140] = 11./7257600.;
		elco[141] = 1./39916800.;
		elco[143] = 4777223./17418240.;
		elco[144] = 1.;
		elco[145] = 83711./55440.;
		elco[146] = 190553./151200.;
		elco[147] = 341747./518400.;
		elco[148] = 139381./604800.;
		elco[149] = 242537./4354560.;
		elco[150] = 1903./201600.;
		elco[151] = 10831./9676800.;
		elco[152] = 11./120960.;
		elco[153] = 1./207360.;
		elco[154] = 1./6652800.;
		elco[155] = 1./479001600.;

		tesco[0] = 0.;
		tesco[1] = 2.;
		tesco[2] = 12.;
		tesco[3] = 1.;
		tesco[4] = 12.;
		tesco[5] = 24.;
		tesco[6] = 2.;
		tesco[7] = 24.;
		tesco[8] = 720./19.;
		tesco[9] = 1.;
		tesco[10] = 720./19.;
		tesco[11] = 160./3.;
		tesco[12] = 6./19.;
		tesco[13] = 160./3.;
		tesco[14] = 60480./863.;
		tesco[15] = 2./27.;
		tesco[16] = 60480./863.;
		tesco[17] = 24192./275.;
		tesco[18] = 12./863.;
		tesco[19] = 24192./275.;
		tesco[20] = 3628800./33953.;
		tesco[21] = 3./1375.;
		tesco[22] = 3628800./33953.;
		tesco[23] = 1036800./8183.;
		tesco[24] = 10./33953.;
		tesco[25] = 1036800./8183.;
		tesco[26] = 479001600./3250433.;
		tesco[27] = 2./57281.;
		tesco[28] = 479001600./3250433.;
		tesco[29] = 788480./4671.;
		tesco[30] = 12./3250433.;
		tesco[31] = 788480./4671.;
		tesco[32] = 2395008.*479001600./6007600823071.;
		tesco[33] = 2./5675265.;
		tesco[34] = 2395008.*479001600./6007600823071.;
		tesco[35] = 0.;

		return;
	}

	static void SetBDFCoefficients(double * const elco, double * const tesco)
	// meth == 2
	{
		elco[0] = 1.;
		elco[1] = 1.;
		elco[13] = 2./3.;
		elco[14] = 1.;
		elco[15] = 1./3.;
		elco[26] = 6./11.;
		elco[27] = 1.;
		elco[28] = 6./11.;
		elco[29] = 1./11.;
		elco[39] = 0.48;
		elco[40] = 1.;
		elco[41] = 7./10.;
		elco[42] = 1./5.;
		elco[43] = 1./50.;
		elco[52] = 60./137.;
		elco[53] = 1.;
		elco[54] = 225./274.;
		elco[55] = 85./274.;
		elco[56] = 15./274.;
		elco[57] = 1./274.;

		tesco[0] = 1.;
		tesco[1] = 2.;
		tesco[2] = 3.;
		tesco[3] = 1.;
		tesco[4] = 9./2.;
		tesco[5] = 6.;
		tesco[6] = 1./2.;
		tesco[7] = 22./3.;
		tesco[8] = 55./6.;
		tesco[9] = 1./6.;
		tesco[10] = 125./12.;
		tesco[11] = 12.5;
		tesco[12] = 1./24.;
		tesco[13] = 137./10.;
		tesco[14] = 959./60.;

		return;
	}

#endif
	protected:
	struct common1 {
		double	ccmax,	// maximum relative change in h__ * EL0 before PJAC is called.
			el0,	//
			h,	// next step size
			hmin,	//
			hmxi,	//
			hu,	// last step size
			rc,	//
			tn;	// last value of t reached

		long	icf,
			ierpj,
			iersl,
			jcur,
			jstart,
			kflag,
			l,
			lewt,
			lacor,
			lsavf,
			lwm,
			msbp,
			mxncf,
			meth,		// current method. 1 means Adams method (nonstiff), 2 means BDF method (stiff)
			miter,		// indicates the corrector iteration method
			maxord,		// maximum order of integration method to be allowed
			maxcor,		// the maximum number of corrector iterations allowed
			n,		// number of equations in current calculations
			nq,		// next method order
			nst,		// number of steps taken so far
			nqu;		// the method order last used (successfully)

		common1()
		{
			ccmax = 0.0;
			el0 = 0.0;
			h = 0.0;
			hmin = 0.0;
			hmxi = 0.0;
			hu = 0.0;
			rc = 0.0;
			tn = 0.0;

			icf = 0;
			ierpj = 0;
			iersl = 0;
			jcur = 0;
			jstart = 0;
			kflag = 0;
			l = 0;
			lewt = 0;
			lacor = 0;
			lsavf = 0;
			lwm = 0;
			meth = 0;
			miter = 0;
			maxord = 0;
			maxcor = 0;
			msbp = 0;
			mxncf = 0;
			n = 0;
			nq = 0;
			nst = 0;
			nqu = 0;
		}
	} dls001;

	struct common2 {
		double pdnorm;
		long jtyp, mused, mxordn, mxords;

		common2()
		{
			pdnorm = 0.0;
			jtyp = 2;
			mused = 0;
			mxordn = 0;
			mxords = 0;
		}
	} dlsa01;

	ErrorMessenger em;
	// TODO Loose ends go here:
	long nyh, liw, lrw, itol, *iwork;
	long itask, istate, iopt;
	double *rwork, *atol, *rtol;
	// ... that's it for loose ends.

	virtual void dsolsy_(double * const x);

	virtual void dprja_(void);
	void dlsoda_(double tout);

	public:
	LivermoreSolver(ODESystem& ODEsys) : ODESolver(ODEsys), iwork(0), rwork(0), itol(1), itask(1), istate(1), iopt(0)
	{
		nyh = ODEsys.nVariables();
		liw = 20 + nyh;
		lrw = 22 + nyh * max((long)16, nyh + 9);
		iwork = new long[liw];
		rwork = new double[lrw];

		for (long i = 0; i < liw; ++i) iwork[i] = 0;
		for (long i = 0; i < lrw; ++i) rwork[i] = 0.;

		rtol = new double(1e-9);
		atol = new double(numeric_limits<double>::epsilon() * 1e3);
	}

	~LivermoreSolver()
	{
		delete [] iwork;
		delete [] rwork;
		delete [] rtol;
		delete [] atol;
	}

	virtual int Solve(const double &);

	void Step(void);

	int InterpolateVariables(const double &);
	int InterpolateDerivatives(const double &, const int, double * const);
};

#endif
