//
// Copyright (c) 2003 Joao Abecasis
//

#include "LivermoreSolver.h"

#include "LinearAlgebra.h"
#include "LinearAlgebra2.h"

using namespace std;
using namespace LSODA::LinearAlgebra;

int LivermoreSolver::InterpolateDerivatives(const double &t, const int order, double * const dky)
{
	long i, j, ic;
	double s;

	/*
	 *  Interpolate solution derivatives.
	 *
	 *  Computes interpolated values of the K-th derivative of the
	 *  dependent variable vector Y, and stores it in DKY.  This routine
	 *  is called within the package with K = 0 and T = TOUT, but may
	 *  also be called by the user for any K up to the current order.
	 *  (See detailed instructions in the usage documentation.)
	 * 
	 *  The computed values in DKY are gotten by interpolation using the
	 *  Nordsieck history array YH.  This array corresponds uniquely to a
	 *  vector-valued polynomial of degree NQCUR or less, and DKY is set
	 *  to the K-th derivative of this polynomial at T.
	 *  The formula for DKY is:
	 *               q
	 *   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
	 *              j=K
	 *
	 *  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
	 *  The above sum is done in reverse order.
	 *
	 *  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are
	 *  communicated by COMMON.
	 *
	 *  Return value is negative if either T or K is out of bounds.
	 *
	 *  The input parameters are:
	 *  T         = value of independent variable where answers are desired
	 *              (normally the same as the T last returned by DLSODA).
	 *              For valid results, T must lie between TCUR - HU and TCUR.
	 *              (See optional outputs for TCUR and HU.)
	 *  K         = integer order of the derivative desired.  K must satisfy
	 *              0 .le. K .le. NQCUR, where NQCUR is the current order
	 *              (see optional outputs).  The capability corresponding
	 *              to K = 0, i.e. computing Y(T), is already provided
	 *              by DLSODA directly.  Since NQCUR .ge. 1, the first
	 *              derivative dy/dt is always available with DINTDY.
	 *  The output parameters are:
	 *  DKY       = a real array of length NEQ containing the computed value
	 *              of the K-th derivative of Y(t).
	 *  Return Value	returned as 0 if K and T were legal,
	 *              -1 if K was illegal, and -2 if T was illegal.
	 */

	if (order < 0 || order > dls001.nq)
	{
		return -1;
	}

	s = dls001.tn - dls001.hu - std::numeric_limits<double>::epsilon() * 100. * (dls001.tn + dls001.hu);
	if ((t - s) * (t - dls001.tn) > 0.) // t must be within hu of tn (error margin 100*epsilon)
	{
		return -2;
	}

	s = (t - dls001.tn) / dls001.h;
	ic = 1;
	if (order != 0)
	{
		for (i = dls001.l - order; i <= dls001.nq; ++i)
		{
			ic *= i;
		}
	}
	for (i = 0; i < dls001.n; ++i)
	{
		dky[i] = ic * rwork[20 + i + (dls001.l - 1) * nyh];
	}
	if (order != dls001.nq)
	{
		for (j = dls001.nq - 1; j >= order; --j)
		{
			ic = 1;
			if (order != 0)
			{
				for (i = j + 1 - order; i <= j; ++i)
				{
					ic *= i;
				}
			}
			for (i = 0; i < dls001.n; ++i)
			{
				dky[i] = ic * rwork[20 + i + j * nyh] + s * dky[i];
			}
		}
		if (order == 0)
		{
			return 0;
		}
	}
	s = pow(dls001.h, -order);
	for (i = 0; i < dls001.n; ++i)
	{
		dky[i] *= s;
	}
	return 0;
}

int LivermoreSolver::InterpolateVariables(const double &t)
{
	// From InterpolateDerivatives, with k = 0 and output to ODEs;
	long i, j;
	double s, *dky;

	dky = new double[dls001.n];

	s = dls001.tn - dls001.hu - std::numeric_limits<double>::epsilon() * 100. * (dls001.tn + dls001.hu);
	if ((t - s) * (t - dls001.tn) > 0.) // t must be within hu of tn (error margin 100*epsilon)
	{
		return -2;
	}

	s = (t - dls001.tn) / dls001.h;
	for (i = 0; i < dls001.n; ++i)
	{
		dky[i] = rwork[20 + i + (dls001.l - 1) * nyh];
	}
	for (j = dls001.nq - 1; j >= 0; --j)
	{
		for (i = 0; i < dls001.n; ++i)
		{
			dky[i] = rwork[20 + i + j * nyh] + s * dky[i];
		}
	}

	for (i = 0; i < dls001.n; ++i)
	{
		ODEs.SetVariable(i, dky[i]);
	}

	delete [] dky;

	return 0;
}

void LivermoreSolver::Step(void)
{
	const double sm1[12] = { .5,.575,.55,.45,.35,.25,.2,.15,.1,.075,.05,.025 };

	long i, j, m, lm1, lm2, ncf, nqm1, nqm2, iret, newq, iredo;
	double r, rh, rm, dm1, dm2, rh1, rh2, del, ddn, pdh, dsm, dup, exm1, exm2, dcon, delp, exdn, rhdn, told, rhsm, exsm, rhup, rate, exup, rh1it, alpha, pnorm;

	static double conit, crate, el[13], elco[156], hold, rmax, tesco[36], cm1[12], cm2[5], pdest, pdlast, ratio;
	static long ialth, ipup, lmax, nqnyh, nslp, icount, irflag;

	/*
	 * Performs one step of the integration of an initial value 
	 * problem for a system of ordinary differential equations.
	 *
	 * Note: DSTODA is independent of the value of the iteration method
	 * indicator MITER, when this is .ne. 0, and hence is independent
	 * of the type of chord method used, or the Jacobian structure.
	 *
	 * Communication with DSTODA is done with the following variables:
	 *
	 * Y      = an array of length .ge. N used as the Y argument in
	 *          all calls to F and JAC.
	 * NEQ    = integer array containing problem size in NEQ(1), and
	 *          passed as the NEQ argument in all calls to F and JAC.
	 * YH     = an NYH by LMAX array containing the dependent variables
	 *          and their approximate scaled derivatives, where
	 *          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
	 *          j-th derivative of Y(i), scaled by H**j/factorial(j)
	 *          (j = 0,1,...,NQ).  On entry for the first step, the first
	 *          two columns of YH must be set from the initial values.
	 * NYH    = a constant integer .ge. N, the first dimension of YH.
	 * YH1    = a one-dimensional array occupying the same space as YH.
	 * EWT    = an array of length N containing multiplicative weights
	 *          for local error measurements.  Local errors in Y(i) are
	 *          compared to 1.0/EWT(i) in various error tests.
	 * SAVF   = an array of working storage, of length N.
	 * ACOR   = a work array of length N, used for the accumulated
	 *          corrections.  On a successful return, ACOR(i) contains
	 *          the estimated one-step local error in Y(i).
	 * WM,IWM = real and integer work arrays associated with matrix
	 *          operations in chord iteration (MITER .ne. 0).
	 * PJAC   = name of routine to evaluate and preprocess Jacobian matrix
	 *          and P = I - H*EL0*Jac, if a chord method is being used.
	 *          It also returns an estimate of norm(Jac) in PDNORM.
	 * SLVS   = name of routine to solve linear system in chord iteration.
	 * CCMAX  = maximum relative change in H*EL0 before PJAC is called.
	 * H      = the step size to be attempted on the next step.
	 *          H is altered by the error control algorithm during the
	 *          problem.  H can be either positive or negative, but its
	 *          sign must remain constant throughout the problem.
	 * HMIN   = the minimum absolute value of the step size H to be used.
	 * HMXI   = inverse of the maximum absolute value of H to be used.
	 *          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
	 *          HMIN and HMXI may be changed at any time, but will not
	 *          take effect until the next change of H is considered.
	 * TN     = the independent variable. TN is updated on each step taken.
	 * JSTART = an integer used for input only, with the following
	 *          values and meanings:
	 *               0  perform the first step.
	 *           .gt.0  take a new step continuing from the last.
	 *              -1  take the next step with a new value of H,
	 *                    N, METH, MITER, and/or matrix parameters.
	 *              -2  take the next step with a new value of H,
	 *                    but with other inputs unchanged.
	 *          On return, JSTART is set to 1 to facilitate continuation.
	 * KFLAG  = a completion code with the following meanings:
	 *               0  the step was succesful.
	 *              -1  the requested error could not be achieved.
	 *              -2  corrector convergence could not be achieved.
	 *              -3  fatal error in PJAC or SLVS.
	 *          A return with KFLAG = -1 or -2 means either
	 *          ABS(H) = HMIN or 10 consecutive failures occurred.
	 *          On a return with KFLAG negative, the values of TN and
	 *          the YH array are as of the beginning of the last
	 *          step, and H is the last step size attempted.
	 * MAXORD = the maximum order of integration method to be allowed.
	 * MAXCOR = the maximum number of corrector iterations allowed.
	 * MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
	 * MXNCF  = maximum number of convergence failures allowed.
	 * METH   = current method.
	 *          METH = 1 means Adams method (nonstiff)
	 *          METH = 2 means BDF method (stiff)
	 *          METH may be reset by DSTODA.
	 * MITER  = corrector iteration method.
	 *          MITER = 0 means functional iteration.
	 *          MITER = JT .gt. 0 means a chord iteration corresponding
	 *          to Jacobian type JT.  (The DLSODA/DLSODAR argument JT is
	 *          communicated here as JTYP, but is not used in DSTODA
	 *          except to load MITER following a method switch.)
	 *          MITER may be reset by DSTODA.
	 * N      = the number of first-order differential equations.
	 */

	dls001.kflag = 0;
	told = dls001.tn;
	ncf = 0;
	dls001.ierpj = 0;
	dls001.iersl = 0;
	dls001.jcur = 0;
	dls001.icf = 0;
	delp = 0.;
	if (dls001.jstart > 0)
	{
		goto L200;
	}
	if (dls001.jstart == -1)
	{
		goto L100;
	}
	if (dls001.jstart == -2)
	{
		goto L160;
	}
	lmax = dls001.maxord + 1;
	dls001.nq = 1;
	dls001.l = 2;
	ialth = 2;
	rmax = 1e4;
	dls001.rc = 0.;
	dls001.el0 = 1.;
	crate = .7;
	hold = dls001.h;
	nslp = 0;
	ipup = dls001.miter;
	iret = 3;
	icount = 20;
	irflag = 0;
	pdest = 0.;
	pdlast = 0.;
	ratio = 5.;
	SetBDFCoefficients(elco, tesco);
	for (i = 0; i < 5; ++i)
	{
		cm2[i] = tesco[i * 3 + 1] * elco[i * 14 + 1];
	}
	SetAdamsCoefficients(elco, tesco);
	for (i = 0; i < 12; ++i)
	{
		cm1[i] = tesco[i * 3 + 1] * elco[i * 14 + 1];
	}
	goto L150;
	L100:
	ipup = dls001.miter;
	lmax = dls001.maxord + 1;
	if (ialth == 1)
	{
		ialth = 2;
	}
	if (dls001.meth == dlsa01.mused)
	{
		goto L160;
	}
	switch (dls001.meth)
	{
		case 1:
			SetAdamsCoefficients(elco, tesco);
			break;
		case 2:
			SetBDFCoefficients(elco, tesco);
			break;
	}
	ialth = dls001.l;
	iret = 1;
	L150:
	for (i = 0; i < dls001.l; ++i)
	{
		el[i] = elco[i + (dls001.nq - 1) * 13];
	}
	nqnyh = dls001.nq * nyh;
	dls001.rc = dls001.rc * el[0] / dls001.el0;
	dls001.el0 = el[0];
	conit = .5 / (dls001.nq + 2);
	switch (iret)
	{
		case 1:  goto L160;
		case 2:  goto L170;
		case 3:  goto L200;
	}
	L160:
	if (dls001.h == hold)
	{
		goto L200;
	}
	rh = dls001.h / hold;
	dls001.h = hold;
	iredo = 3;
	goto L175;
	L170:
	rh = max(rh, dls001.hmin / fabs(dls001.h));
	L175:
	rh = min(rh,rmax);
	rh /= max(1., fabs(dls001.h) * dls001.hmxi * rh);
	if (dls001.meth == 2)
	{
		goto L178;
	}
	irflag = 0;
	pdh = max(fabs(dls001.h) * pdlast,1e-6);
	if (rh * pdh * 1.00001 < sm1[dls001.nq - 1])
	{
		goto L178;
	}
	rh = sm1[dls001.nq - 1] / pdh;
	irflag = 1;
	L178:
	r = 1.;
	for (j = 1; j < dls001.l; ++j)
	{
		r *= rh;
		for (i = 0; i < dls001.n; ++i)
		{
			rwork[20 + i + j * nyh] *= r;
		}
	}
	dls001.h *= rh;
	dls001.rc *= rh;
	ialth = dls001.l;
	if (iredo == 0)
	{
		goto L690;
	}
	L200:
	if (fabs(dls001.rc - 1.) > dls001.ccmax)
	{
		ipup = dls001.miter;
	}
	if (dls001.nst >= nslp + dls001.msbp)
	{
		ipup = dls001.miter;
	}
	dls001.tn += dls001.h;
	for (j = nqnyh - nyh; j >= 0; j -= nyh)
	{
		for (i = j; i < nqnyh; ++i)
		{
			rwork[20 + i] += rwork[20 + i + nyh];
		}
	}
	pnorm = VectorMaxNorm(dls001.n, &rwork[20], &rwork[dls001.lewt - 1]);
	L220:
	m = 0;
	rate = 0.;
	del = 0.;
	for (i = 0; i < dls001.n; ++i)
	{
		ODEs.SetVariable(i, rwork[20 + i]);
	}
	ODEs.SetTime(dls001.tn);
	for (i = 0; i < dls001.n; ++i)
	{
		rwork[dls001.lsavf + i - 1] = ODEs.GetDerivative(i);
	}
	if (ipup > 0)
	{
		this->dprja_();
		ipup = 0;
		dls001.rc = 1.;
		nslp = dls001.nst;
		crate = .7;
		if (dls001.ierpj != 0)
		{
			goto L430;
		}
	}
	for (i = 0; i < dls001.n; ++i)
	{
		rwork[dls001.lacor + i - 1] = 0.;
	}
	L270:
	if (dls001.miter != 0)
	{
		goto L350;
	}
	del = 0.;
	for (i = 0; i < dls001.n; ++i)
	{
		rwork[dls001.lsavf + i - 1] = dls001.h * rwork[dls001.lsavf + i - 1] - rwork[20 + i + nyh];
		del = max(del, fabs(rwork[dls001.lsavf + i - 1] - rwork[dls001.lacor + i - 1]) * rwork[dls001.lewt + i - 1]);
	}
	for (i = 0; i < dls001.n; ++i)
	{
		ODEs.SetVariable(i, rwork[20 + i] + el[0] * rwork[dls001.lsavf + i - 1]);
		rwork[dls001.lacor + i - 1] = rwork[dls001.lsavf + i - 1];
	}
	goto L400;
	L350:
	for (i = 0; i < dls001.n; ++i)
	{
		ODEs.SetVariable(i, dls001.h * rwork[dls001.lsavf + i - 1] - (rwork[20 + i + nyh] + rwork[dls001.lacor + i - 1]));
	}
	this->dsolsy_(&rwork[dls001.lsavf - 1]);
	if (dls001.iersl < 0)
	{
		goto L430;
	}
	if (dls001.iersl > 0)
	{
		goto L410;
	}
	del = 0.;
	for (i = 0; i < dls001.n; ++i)
	{
		del = max(del, fabs(ODEs.GetVariable(i)) * rwork[dls001.lewt + i - 1]);
	}
	for (i = 0; i < dls001.n; ++i)
	{
		rwork[dls001.lacor + i - 1] += ODEs.GetVariable(i);
		ODEs.SetVariable(i, rwork[20 + i] + el[0] * rwork[dls001.lacor + i - 1]);
	}
	L400:
	if (del <= pnorm * 100. * numeric_limits<double>::epsilon())
	{
		goto L450;
	}
	if (m == 0 && dls001.meth == 1)
	{
		goto L405;
	}
	if (m == 0)
	{
		goto L402;
	}
	rm = 1024.;
	if (del <= delp * 1024.)
	{
		rm = del / delp;
	}
	rate = max(rate,rm);
	crate = max(crate * .2,rm);
	L402:
	dcon = del * min(1., crate * 1.5) / (tesco[dls001.nq * 3 - 2] * conit);
	if (dcon > 1.)
	{
		goto L405;
	}
	pdest = max(pdest, rate / fabs(dls001.h * el[0]));
	if (pdest != 0.)
	{
		pdlast = pdest;
	}
	goto L450;
	L405:
	++m;
	if (m == dls001.maxcor)
	{
		goto L410;
	}
	if (m >= 2 && del > delp * 2.)
	{
		goto L410;
	}
	delp = del;
	ODEs.SetTime(dls001.tn);
	for (i = 0; i < dls001.n; ++i)
	{
		rwork[dls001.lsavf + i - 1] = ODEs.GetDerivative(i);
	}
	goto L270;
	L410:
	if (dls001.miter == 0 || dls001.jcur == 1)
	{
		goto L430;
	}
	dls001.icf = 1;
	ipup = dls001.miter;
	goto L220;
	L430:
	dls001.icf = 2;
	++ncf;
	rmax = 2.;
	dls001.tn = told;
	for (j = nqnyh - nyh; j >= 0; j -= nyh)
	{
		for (i = j; i < nqnyh; ++i)
		{
			rwork[20 + i] -= rwork[20 + i + nyh];
		}
	}
	if (dls001.ierpj < 0 || dls001.iersl < 0)
	{
		goto L680;
	}
	if (fabs(dls001.h) <= dls001.hmin * 1.00001)
	{
		goto L670;
	}
	if (ncf == dls001.mxncf)
	{
		goto L670;
	}
	rh = .25;
	ipup = dls001.miter;
	iredo = 1;
	goto L170;
	L450:
	dls001.jcur = 0;
	if (m == 0)
	{
		dsm = del / tesco[dls001.nq * 3 - 2];
	}
	if (m > 0)
	{
		dsm = VectorMaxNorm(dls001.n, &rwork[dls001.lacor - 1], &rwork[dls001.lewt - 1]) / tesco[dls001.nq * 3 - 2];
	}
	if (dsm > 1.)
	{
		goto L500;
	}
	dls001.kflag = 0;
	iredo = 0;
	++dls001.nst;
	dls001.hu = dls001.h;
	dls001.nqu = dls001.nq;
	dlsa01.mused = dls001.meth;
	for (j = 0; j < dls001.l; ++j)
	{
		for (i = 0; i < dls001.n; ++i)
		{
			rwork[20 + i + j * nyh] += el[j] * rwork[dls001.lacor + i - 1];
		}
	}
	--icount;
	if (icount >= 0)
	{
		goto L488;
	}
	if (dls001.meth == 2)
	{
		goto L480;
	}
	if (dls001.nq > 5)
	{
		goto L488;
	}
	if (dsm > pnorm * 100. * numeric_limits<double>::epsilon() && pdest != 0.)
	{
		goto L470;
	}
	if (irflag == 0)
	{
		goto L488;
	}
	rh2 = 2.;
	nqm2 = min(dls001.nq,dlsa01.mxords);
	goto L478;
	L470:
	exsm = 1. / dls001.l;
	rh1 = 1. / (pow(dsm, exsm) * 1.2 + 1.2e-6);
	rh1it = rh1 * 2.;
	pdh = pdlast * fabs(dls001.h);
	if (pdh * rh1 > 1e-5)
	{
		rh1it = sm1[dls001.nq - 1] / pdh;
	}
	rh1 = min(rh1,rh1it);
	if (dls001.nq <= dlsa01.mxords)
	{
		goto L474;
	}
	nqm2 = dlsa01.mxords;
	lm2 = dlsa01.mxords + 1;
	exm2 = 1. / lm2;
	dm2 = VectorMaxNorm(dls001.n, &rwork[20 + lm2 * nyh], &rwork[dls001.lewt - 1]) / cm2[dlsa01.mxords - 1];
	rh2 = 1. / (pow(dm2, exm2) * 1.2 + 1.2e-6);
	goto L476;
	L474:
	dm2 = dsm * (cm1[dls001.nq - 1] / cm2[dls001.nq - 1]);
	rh2 = 1. / (pow(dm2, exsm) * 1.2 + 1.2e-6);
	nqm2 = dls001.nq;
	L476:
	if (rh2 < ratio * rh1)
	{
		goto L488;
	}
	L478:
	rh = rh2;
	icount = 20;
	dls001.meth = 2;
	dls001.miter = dlsa01.jtyp;
	pdlast = 0.;
	dls001.nq = nqm2;
	dls001.l = dls001.nq + 1;
	goto L170;
	L480:
	exsm = 1. / dls001.l;
	if (dlsa01.mxordn >= dls001.nq)
	{
		goto L484;
	}
	nqm1 = dlsa01.mxordn;
	lm1 = dlsa01.mxordn + 1;
	exm1 = 1. / lm1;
	dm1 = VectorMaxNorm(dls001.n, &rwork[20 + lm1 * nyh], &rwork[dls001.lewt - 1]) / cm1[
		dlsa01.mxordn - 1];
	rh1 = 1. / (pow(dm1, exm1) * 1.2 + 1.2e-6);
	goto L486;
	L484:
	dm1 = dsm * (cm2[dls001.nq - 1] / cm1[dls001.nq - 1]);
	rh1 = 1. / (pow(dm1, exsm) * 1.2 + 1.2e-6);
	nqm1 = dls001.nq;
	exm1 = exsm;
	L486:
	rh1it = rh1 * 2.;
	pdh = dlsa01.pdnorm * fabs(dls001.h);
	if (pdh * rh1 > 1e-5)
	{
		rh1it = sm1[nqm1 - 1] / pdh;
	}
	rh1 = min(rh1,rh1it);
	rh2 = 1. / (pow(dsm, exsm) * 1.2 + 1.2e-6);
	if (rh1 * ratio < rh2 * 5.)
	{
		goto L488;
	}
	alpha = max(.001,rh1);
	dm1 = pow(alpha, exm1) * dm1;
	if (dm1 <= numeric_limits<double>::epsilon() * 1e3 * pnorm)
	{
		goto L488;
	}
	rh = rh1;
	icount = 20;
	dls001.meth = 1;
	dls001.miter = 0;
	pdlast = 0.;
	dls001.nq = nqm1;
	dls001.l = dls001.nq + 1;
	goto L170;
	L488:
	--ialth;
	if (ialth == 0)
	{
		goto L520;
	}
	if (ialth > 1)
	{
		goto L700;
	}
	if (dls001.l == lmax)
	{
		goto L700;
	}
	for (i = 0; i < dls001.n; ++i)
	{
		rwork[20 + i + (lmax - 1) * nyh] = rwork[dls001.lacor + i - 1];
	}
	goto L700;
	L500:
	--dls001.kflag;
	dls001.tn = told;
	for (j = nqnyh - nyh; j >= 0; j -= nyh)
	{
		for (i = j; i < nqnyh; ++i)
		{
			rwork[20 + i] -= rwork[20 + i + nyh];
		}
	}
	rmax = 2.;
	if (fabs(dls001.h) <= dls001.hmin * 1.00001)
	{
		goto L660;
	}
	if (dls001.kflag <= -3)
	{
		goto L640;
	}
	iredo = 2;
	rhup = 0.;
	goto L540;
	L520:
	rhup = 0.;
	if (dls001.l == lmax)
	{
		goto L540;
	}
	for (i = 0; i < dls001.n; ++i)
	{
		rwork[dls001.lsavf + i - 1] = rwork[dls001.lacor + i - 1] - rwork[20 + i + (lmax - 1) * nyh];
	}
	dup = VectorMaxNorm(dls001.n, &rwork[dls001.lsavf - 1], &rwork[dls001.lewt - 1]) / tesco[dls001.nq * 3 - 1];
	exup = 1. / (dls001.l + 1);
	rhup = 1. / (pow(dup, exup) * 1.4 + 1.4e-6);
	L540:
	exsm = 1. / dls001.l;
	rhsm = 1. / (pow(dsm, exsm) * 1.2 + 1.2e-6);
	rhdn = 0.;
	if (dls001.nq == 1)
	{
		goto L550;
	}
	ddn = VectorMaxNorm(dls001.n, &rwork[20 + (dls001.l - 1) * nyh], &rwork[dls001.lewt - 1]) /
		tesco[dls001.nq * 3 - 3];
	exdn = 1. / dls001.nq;
	rhdn = 1. / (pow(ddn, exdn) * 1.3 + 1.3e-6);
	L550:
	if (dls001.meth == 2)
	{
		goto L560;
	}
	pdh = max(fabs(dls001.h) * pdlast,1e-6);
	if (dls001.l < lmax)
	{
		rhup = min(rhup, sm1[dls001.l - 1] / pdh);
	}
	;
	rhsm = min(rhsm, sm1[dls001.nq - 1] / pdh);
	if (dls001.nq > 1)
	{
		rhdn = min(rhdn, sm1[dls001.nq - 2] / pdh);
	}
	pdest = 0.;
	L560:
	if (rhsm >= rhup)
	{
		goto L570;
	}
	if (rhup > rhdn)
	{
		goto L590;
	}
	goto L580;
	L570:
	if (rhsm < rhdn)
	{
		goto L580;
	}
	newq = dls001.nq;
	rh = rhsm;
	goto L620;
	L580:
	newq = dls001.nq - 1;
	rh = rhdn;
	if (dls001.kflag < 0 && rh > 1.)
	{
		rh = 1.;
	}
	goto L620;
	L590:
	newq = dls001.l;
	rh = rhup;
	if (rh < 1.1)
	{
		goto L610;
	}
	r = el[dls001.l - 1] / dls001.l;
	for (i = 0; i < dls001.n; ++i)
	{
		rwork[20 + i + (newq + 1) * nyh] = rwork[dls001.lacor + i - 1] * r;
	}
	goto L630;
	L610:
	ialth = 3;
	goto L700;
	L620:
	if (dls001.meth == 2)
	{
		goto L622;
	}
	if (rh * pdh * 1.00001 >= sm1[newq - 1])
	{
		goto L625;
	}
	L622:
	if (dls001.kflag == 0 && rh < 1.1)
	{
		goto L610;
	}
	L625:
	if (dls001.kflag <= -2)
	{
		rh = min(rh,.2);
	}
	if (newq == dls001.nq)
	{
		goto L170;
	}
	L630:
	dls001.nq = newq;
	dls001.l = dls001.nq + 1;
	iret = 2;
	goto L150;
	L640:
	if (dls001.kflag == -10)
	{
		goto L660;
	}
	rh = .1;
	rh = max(dls001.hmin / fabs(dls001.h),rh);
	dls001.h *= rh;
	for (i = 0; i < dls001.n; ++i)
	{
		ODEs.SetVariable(i, rwork[20 + i]);
	}
	ODEs.SetTime(dls001.tn);
	for (i = 0; i < dls001.n; ++i)
	{
		rwork[dls001.lsavf + i - 1] = ODEs.GetDerivative(i);
	}
	for (i = 0; i < dls001.n; ++i)
	{
		rwork[20 + i + nyh] = dls001.h * rwork[dls001.lsavf + i - 1];
	}
	ipup = dls001.miter;
	ialth = 5;
	if (dls001.nq == 1)
	{
		goto L200;
	}
	dls001.nq = 1;
	dls001.l = 2;
	iret = 3;
	goto L150;
	L660:
	dls001.kflag = -1;
	goto L720;
	L670:
	dls001.kflag = -2;
	goto L720;
	L680:
	dls001.kflag = -3;
	goto L720;
	L690:
	rmax = 10.;
	L700:
	r = 1. / tesco[dls001.nqu * 3 - 2];
	for (i = 0; i < dls001.n; ++i)
	{
		rwork[dls001.lacor + i - 1] *= r;
	}
	L720:
	hold = dls001.h;
	dls001.jstart = 1;
	return;
}

void LivermoreSolver::dsolsy_(double * const x)
{
	long i, ml, mu, meband;
	double r, di, hl0, phl0;

	/*  Linear system solver.
	 *
	 *  This routine manages the solution of the linear system arising from
	 *  a chord iteration.  It is called if MITER .ne. 0.
	 *  If MITER is 1 or 2, it calls DGESL to accomplish this.
	 *  If MITER = 3 it updates the coefficient h*EL0 in the diagonal
	 *  matrix, and then computes the solution.
	 *  If MITER is 4 or 5, it calls DGBSL.
	 *  Communication with DSOLSY uses the following variables:
	 *  WM    = real work space containing the inverse diagonal matrix if
	 *          MITER = 3 and the LU decomposition of the matrix otherwise.
	 *          Storage of matrix elements starts at WM(3).
	 *          WM also contains the following matrix-related data:
	 *          WM(1) = SQRT(UROUND) (not used here),
	 *          WM(2) = HL0, the previous value of h*EL0, used if MITER = 3.
	 *  IWM   = integer work space containing pivot information, starting at
	 *          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
	 *          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
	 *  X     = the right-hand side vector on input, and the solution vector
	 *          on output, of length N.
	 *  IERSL = output flag (in COMMON).  IERSL = 0 if no trouble occurred.
	 *          IERSL = 1 if a singular matrix arose with MITER = 3.
	 *
	 *  This routine also uses the COMMON variables EL0, H, MITER, and N.
	 *
	 *  Routines called: dgbsl_, dgesl_
	 */

	dls001.iersl = 0;
	switch (dls001.miter)
	{
		case 1:
		case 2:
			dgesl_(rwork + dls001.lwm + 1, dls001.n, dls001.n, &iwork[20], x, 0);
			return;
		case 3:
			phl0 = rwork[dls001.lwm];
			hl0 = dls001.h * dls001.el0;
			rwork[dls001.lwm] = hl0;
			if (hl0 != phl0)
			{
				r = hl0 / phl0;
				for (i = 0; i < dls001.n; ++i)
				{
					di = 1. - r * (1. - 1. / rwork[dls001.lwm + i + 1]);
					if (fabs(di) == 0.)
					{
						dls001.iersl = 1;
						return;
					}
					rwork[dls001.lwm + i + 1] = 1. / di;
				}
			}
			for (i = 0; i < dls001.n; ++i)
			{
				x[i] = rwork[dls001.lwm + i + 1] * x[i];
			}
			return;
		case 4:
		case 5:
			ml = iwork[0];
			mu = iwork[1];
			meband = (ml << 1) + mu + 1;
			dgbsl_(&rwork[dls001.lwm + 1], meband, dls001.n, &ml, &mu, &iwork[20], x, 0);
			return;
	}
}

void LivermoreSolver::dprja_(void)
{
	long i__, j, i1, i2, j1, ii, jj, ml, mu, np1, mba, ier, meb1, lenp, mband, meband;
	double r__, r0, yi, yj, hl0, fac, con, yjj, srur;

	/*
	 * Called by LivermoreSolver::Step() to compute and process the matrix
	 * P = I - H*EL(1)*J , where J is an approximation to the Jacobian.
	 * Here J is computed by the user-supplied routine JAC if
	 * MITER = 1 or 4 or by finite differencing if MITER = 2 or 5.
	 * J, scaled by -H*EL(1), is stored in WM.  Then the norm of J (the
	 * matrix norm consistent with the weighted max-norm on vectors given
	 * by DMNORM) is computed, and J is overwritten by P.  P is then
	 * subjected to LU decomposition in preparation for later solution
	 * of linear systems with P as coefficient matrix.  This is done
	 * by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
	 * 
	 * In addition to variables described previously, communication
	 * with DPRJA uses the following:
	 * Y     = array containing predicted values on entry.
	 * FTEM  = work array of length N (ACOR in DSTODA).
	 * SAVF  = array containing f evaluated at predicted Y.
	 * WM    = real work space for matrices.  On output it contains the
	 *         LU decomposition of P.
	 *         Storage of matrix elements starts at WM(3).
	 *         WM also contains the following matrix-related data:
	 *         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
	 * IWM   = integer work space containing pivot information, starting at
	 *         IWM(21).   IWM also contains the band parameters
	 *         ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
	 * EL0   = EL(1) (input).
	 * PDNORM= norm of Jacobian matrix. (Output).
	 * IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
	 *         P matrix found to be singular.
	 * JCUR  = output flag = 1 to indicate that the Jacobian matrix
	 *         (or approximation) is now current.
	 * This routine also uses the Common variables EL0, H, TN, UROUND,
	 * MITER, N, NFE, and NJE.
	 */

	dls001.ierpj = 0;
	dls001.jcur = 1;
	hl0 = dls001.h * dls001.el0;
	switch (dls001.miter)
	{
		case 1:  goto L100;
		case 2:  goto L200;
		case 3:  goto L300;
		case 4:  goto L400;
		case 5:  goto L500;
	}
	/* If MITER = 1, call JAC and multiply by scalar. ----------------------- */
	L100:
	lenp = dls001.n * dls001.n;
	for (i__ = 1; i__ <= lenp; ++i__)
	{
		rwork[dls001.lwm + i__] = 0.;
	}
	// TODO: CALL JACOBIAN (TIME = TN, Y, 0, 0, WM(N,*), N)
	ODEs.Jacobian(&rwork[dls001.lwm + 2]);
	con = -hl0;
	for (i__ = 1; i__ <= lenp; ++i__)
	{
		rwork[dls001.lwm + i__ ] *= con;
	}
	goto L240;
	/* If MITER = 2, make N calls to F to approximate J. -------------------- */
	L200:
	fac = VectorMaxNorm(dls001.n, &rwork[dls001.lsavf - 1], &rwork[dls001.lewt - 1]);
	r0 = fabs(dls001.h) * 1e3 * std::numeric_limits<double>::epsilon() * dls001.n * fac;
	if (r0 == 0.)
	{
		r0 = 1.;
	}
	srur = rwork[dls001.lwm-1];
	j1 = 2;
	ODEs.SetTime(dls001.tn);
	for (j = 1; j <= dls001.n; ++j)
	{
		yj = ODEs.GetVariable(j - 1);
		/* Computing MAX */
		r__ = max(srur * fabs(yj), r0 / rwork[dls001.lewt + j - 2]);
		ODEs.SetVariable(j - 1, yj + r__);
		fac = -hl0 / r__;
		for (i__ = 1; i__ <= dls001.n; ++i__)
		{
			rwork[dls001.lwm + i__ + j1 - 2] = (ODEs.GetDerivative(i__ - 1) - rwork[dls001.lsavf + i__ - 2]) * fac;
		}
		ODEs.SetVariable(j - 1, yj);
		j1 += dls001.n;
	}
	L240:
	/* Compute norm of Jacobian. -------------------------------------------- */
	dlsa01.pdnorm = MatrixMaxNorm(dls001.n, &rwork[dls001.lwm + 1], &rwork[dls001.lewt - 1]) / fabs(hl0);
	/* Add identity matrix. ------------------------------------------------- */
	j = 3;
	np1 = dls001.n + 1;
	for (i__ = 1; i__ <= dls001.n; ++i__)
	{
		rwork[dls001.lwm + j -2] += 1.;
		/* L250: */
		j += np1;
	}
	/* Do LU decomposition on P. -------------------------------------------- */
	dgefa_(&rwork[dls001.lwm + 1], dls001.n, dls001.n, &iwork[20], &ier);
	if (ier != 0)
	{
		dls001.ierpj = 1;
	}
	return;
	/* Dummy block only, since MITER is never 3 in this routine. ------------ */
	L300:
	return;
	/* If MITER = 4, call JAC and multiply by scalar. ----------------------- */
	L400:
	ml = iwork[0];
	mu = iwork[1];
	mband = ml + mu + 1;
	meband = mband + ml;
	lenp = meband * dls001.n;
	for (i__ = 1; i__ <= lenp; ++i__)
	{
		rwork[dls001.lwm + i__] = 0.;
	}
	// TODO: CALL JACOBIAN (TIME = TN, Y, ML, MU, WM(meband,*), meband)
	ODEs.Jacobian(&rwork[dls001.lwm + ml + 2]);
	con = -hl0;
	for (i__ = 1; i__ <= lenp; ++i__)
	{
		rwork[dls001.lwm + i__] *= con;
	}
	goto L570;
	/* If MITER = 5, make MBAND calls to F to approximate J. ---------------- */
	L500:
	ml = iwork[0];
	mu = iwork[1];
	mband = ml + mu + 1;
	mba = min(mband,dls001.n);
	meband = mband + ml;
	meb1 = meband - 1;
	srur = rwork[dls001.l - 1];
	fac = VectorMaxNorm(dls001.n, &rwork[dls001.lsavf - 1], &rwork[dls001.lewt - 1]);
	r0 = fabs(dls001.h) * 1e3 * std::numeric_limits<double>::epsilon() * dls001.n * fac;
	if (r0 == 0.)
	{
		r0 = 1.;
	}
	ODEs.SetTime(dls001.tn);
	for (j = 1; j <= mba; ++j)
	{
		for (i__ = j; mband < 0 ? i__ >= dls001.n : i__ <= dls001.n; i__ += mband)
		{
			yi = ODEs.GetVariable(i__ - 1);
			/* Computing MAX */
			r__ = max(srur * fabs(yi), r0 / rwork[dls001.lewt + i__ - 2]);
			ODEs.SetVariable(i__ - 1, yi + r__);
		}
		for (jj = j; mband < 0 ? jj >= dls001.n : jj <= dls001.n; jj += mband)
		{
			ODEs.SetVariable(jj - 1, rwork[jj + 19]);
			yjj = ODEs.GetVariable(jj - 1);
			/* Computing MAX */
			r__ = max(srur * fabs(yjj), r0 / rwork[dls001.lewt + jj - 2]);
			fac = -hl0 / r__;
			/* Computing MAX */
			i1 = max(jj - mu, (long) 1);
			/* Computing MIN */
			i2 = min(jj + ml,dls001.n);
			ii = jj * meb1 - ml + 2;
			for (i__ = i1; i__ <= i2; ++i__)
			{
				rwork[dls001.lwm + ii + i__ - 2] = (ODEs.GetDerivative(i__ - 1) - rwork[dls001.lsavf + i__-2]) * fac;
			}
		}
	}
	L570:
	/* Compute norm of Jacobian. -------------------------------------------- */
	dlsa01.pdnorm = BandedMatrixMaxNorm(dls001.n, &rwork[dls001.lwm + ml + 1], meband, ml, mu, &rwork[dls001.lewt - 1]) / fabs(hl0);
	/* Add identity matrix. ------------------------------------------------- */
	ii = mband + 2;
	for (i__ = 1; i__ <= dls001.n; ++i__)
	{
		rwork[dls001.lwm + ii - 2] += 1.;
		ii += meband;
	}
	/* Do LU decomposition of P. -------------------------------------------- */
	dgbfa_(&rwork[dls001.lwm + 1], meband, dls001.n, &ml, &mu, &iwork[20], &ier);
	if (ier != 0)
	{
		dls001.ierpj = 1;
	}
	return;
}

void LivermoreSolver::dlsoda_(double tout)
{
	/* Initialized data */
	const long mord[2] = {12, 5};
	const long mxstp0 = 5000;
	const long mxhnl0 = 10;

	/* Local variables */
	long i, i1, i2, ml, mu, lf0, kgo, imxer, len1, len2, len1c, len1n, len1s, leniw, lenwm, lenrw, leniwc, lenrwc;
	double h0, w0, rh, tp, big, ayi, hmx, tol, sum, hmax, ewti, size, atoli, tcrit, tdist, rtoli, tolsf, tnext;
	std::string Message;
	bool ihit;

	static long init(0), insufi, insufr, ixpr, mxhnil, mxstep, nhnil, nslast;
	static double tsw;

	if (istate < 1 || istate > 3)
	{
		goto L601;
	}
	if (itask < 1 || itask > 5)
	{
		goto L602;
	}
	if (istate == 1)
	{
		goto L10;
	}
	if (init == 0)
	{
		goto L603;
	}
	if (istate == 2)
	{
		goto L200;
	}
	goto L20;
	L10:
	init = 0;
	if (tout == ODEs.GetTime())
	{
		return;
	}
	L20:
	if (istate == 1)
	{
		goto L25;
	}
	if (ODEs.nVariables() > dls001.n)
	{
		goto L605;
	}
	L25:
	dls001.n = ODEs.nVariables();
	if (itol < 1 || itol > 4)
	{
		goto L606;
	}
	if (iopt < 0 || iopt > 1)
	{
		goto L607;
	}
	if (dlsa01.jtyp == 3 || dlsa01.jtyp < 1 || dlsa01.jtyp > 5)
	{
		goto L608;
	}
	if (dlsa01.jtyp <= 2)
	{
		goto L30;
	}
	ml = iwork[0];
	mu = iwork[1];
	if (ml < 0 || ml >= dls001.n)
	{
		goto L609;
	}
	if (mu < 0 || mu >= dls001.n)
	{
		goto L610;
	}
	L30:
	if (iopt == 1)
	{
		goto L40;
	}
	ixpr = 0;
	mxstep = mxstp0;
	mxhnil = mxhnl0;
	dls001.hmxi = 0.;
	dls001.hmin = 0.;
	if (istate != 1)
	{
		goto L60;
	}
	h0 = 0.;
	dlsa01.mxordn = mord[0];
	dlsa01.mxords = mord[1];
	goto L60;
	L40:
	ixpr = iwork[4];
	if (ixpr < 0 || ixpr > 1)
	{
		goto L611;
	}
	mxstep = iwork[5];
	if (mxstep < 0)
	{
		goto L612;
	}
	if (mxstep == 0)
	{
		mxstep = mxstp0;
	}
	mxhnil = iwork[6];
	if (mxhnil < 0)
	{
		goto L613;
	}
	if (mxhnil == 0)
	{
		mxhnil = mxhnl0;
	}
	if (istate != 1)
	{
		goto L50;
	}
	h0 = rwork[4];
	dlsa01.mxordn = iwork[7];
	if (dlsa01.mxordn < 0)
	{
		goto L628;
	}
	if (dlsa01.mxordn == 0)
	{
		dlsa01.mxordn = 100;
	}
	dlsa01.mxordn = min(dlsa01.mxordn, mord[0]);
	dlsa01.mxords = iwork[8];
	if (dlsa01.mxords < 0)
	{
		goto L629;
	}
	if (dlsa01.mxords == 0)
	{
		dlsa01.mxords = 100;
	}
	dlsa01.mxords = min(dlsa01.mxords, mord[1]);
	if ((tout - ODEs.GetTime()) * h0 < 0.)
	{
		goto L614;
	}
	L50:
	hmax = rwork[5];
	if (hmax < 0.)
	{
		goto L615;
	}
	dls001.hmxi = 0.;
	if (hmax > 0.)
	{
		dls001.hmxi = 1. / hmax;
	}
	dls001.hmin = rwork[6];
	if (dls001.hmin < 0.)
	{
		goto L616;
	}
	L60:
	if (istate == 1)
	{
		dls001.meth = 1;
	}
	if (istate == 1)
	{
		nyh = dls001.n;
	}
	len1n = (dlsa01.mxordn + 1) * nyh + 20;
	len1s = (dlsa01.mxords + 1) * nyh + 20;
	dls001.lwm = len1s + 1;
	if (dlsa01.jtyp <= 2)
	{
		lenwm = dls001.n * dls001.n + 2;
	}
	if (dlsa01.jtyp >= 4)
	{
		lenwm = ((ml << 1) + mu + 1) * dls001.n + 2;
	}
	len1s += lenwm;
	len1c = len1n;
	if (dls001.meth == 2)
	{
		len1c = len1s;
	}
	len1 = max(len1n, len1s);
	len2 = dls001.n * 3;
	lenrw = len1 + len2;
	lenrwc = len1c + len2;
	iwork[16] = lenrw;
	leniw = dls001.n + 20;
	leniwc = 20;
	if (dls001.meth == 2)
	{
		leniwc = leniw;
	}
	iwork[17] = leniw;
	if (istate == 1 && lrw < lenrwc)
	{
		goto L617;
	}
	if (istate == 1 && liw < leniwc)
	{
		goto L618;
	}
	if (istate == 3 && lrw < lenrwc)
	{
		goto L550;
	}
	if (istate == 3 && liw < leniwc)
	{
		goto L555;
	}
	dls001.lewt = len1 + 1;
	insufr = 0;
	if (lrw >= lenrw)
	{
		goto L65;
	}
	insufr = 2;
	dls001.lewt = len1c + 1;
	Message = "DLSODA - Warning.. RWORK length is sufficient for now, but may not be later. Integration will proceed anyway. Length needed is LENRW = I1, while LRW = I2.";
	em.Message(Message, 103, 0, 2, lenrw, lrw, 0, 0., 0.);
	L65:
	dls001.lsavf = dls001.lewt + dls001.n;
	dls001.lacor = dls001.lsavf + dls001.n;
	insufi = 0;
	if (liw >= leniw)
	{
		goto L70;
	}
	insufi = 2;
	Message = "DLSODA - Warning.. IWORK length is sufficient for now, but may not be later. Integration will proceed anyway. Length needed is LENIW = I1, while LIW = I2.";
	em.Message(Message, 104, 0, 2, leniw, liw, 0, 0., 0.);
	L70:
	rtoli = rtol[0];
	atoli = atol[0];
	for (i = 0; i < dls001.n; ++i)
	{
		if (itol >= 3)
		{
			rtoli = rtol[i];
		}
		if (itol == 2 || itol == 4)
		{
			atoli = atol[i];
		}
		if (rtoli < 0.)
		{
			goto L619;
		}
		if (atoli < 0.)
		{
			goto L620;
		}
	}
	if (istate == 1)
	{
		goto L100;
	}
	dls001.jstart = -1;
	if (dls001.n == nyh)
	{
		goto L200;
	}
	i1 = 20 + dls001.l * nyh;
	i2 = 19 + (dls001.maxord + 1) * nyh;
	if (i1 > i2)
	{
		goto L200;
	}
	for (i = i1; i <= i2; ++i)
	{
		rwork[i] = 0.;
	}
	goto L200;
	L100:
	dls001.tn = ODEs.GetTime();
	tsw = ODEs.GetTime();
	dls001.maxord = dlsa01.mxordn;
	if (itask != 4 && itask != 5)
	{
		goto L110;
	}
	tcrit = rwork[0];
	if ((tcrit - tout) * (tout - ODEs.GetTime()) < 0.)
	{
		goto L625;
	}
	if (h0 != 0. && (ODEs.GetTime() + h0 - tcrit) * h0 > 0.)
	{
		h0 = tcrit - ODEs.GetTime();
	}
	L110:
	dls001.jstart = 0;
	nhnil = 0;
	dls001.nst = 0;
	nslast = 0;
	dls001.hu = 0.;
	dls001.nqu = 0;
	dlsa01.mused = 0;
	dls001.miter = 0;
	dls001.ccmax = .3;
	dls001.maxcor = 3;
	dls001.msbp = 20;
	dls001.mxncf = 10;
	lf0 = 20 + nyh;
	for (i = 0; i < ODEs.nVariables(); ++i)
	{
		rwork[lf0 + i] = ODEs.GetDerivative(i);
	}
	for (i = 0; i < dls001.n; ++i)
	{
		rwork[i + 20] = ODEs.GetVariable(i);
	}
	dls001.nq = 1;
	dls001.h = 1.;
	ODEs.GetErrorWeights(itol, rtol, atol, &rwork[20], &rwork[dls001.lewt-1]);
	for (i = 0; i < dls001.n; ++i)
	{
		if (rwork[i + dls001.lewt - 1] <= 0.)
		{
			goto L621;
		}
		rwork[i + dls001.lewt - 1] = 1. / rwork[i + dls001.lewt - 1];
	}
	if (h0 != 0.)
	{
		goto L180;
	}
	tdist = fabs(tout - ODEs.GetTime());
	/* Computing MAX */
	w0 = max(fabs(ODEs.GetTime()), fabs(tout));
	if (tdist < numeric_limits<double>::epsilon() * 2. * w0)
	{
		goto L622;
	}
	tol = rtol[0];
	if (itol <= 2)
	{
		goto L140;
	}
	for (i = 0; i < dls001.n; ++i)
	{
		/* Computing MAX */
		tol = max(tol, rtol[i]);
	}
	L140:
	if (tol > 0.)
	{
		goto L160;
	}
	atoli = atol[0];
	for (i = 0; i < dls001.n; ++i)
	{
		if (itol == 2 || itol == 4)
		{
			atoli = atol[i];
		}
		ayi = fabs(ODEs.GetVariable(i));
		if (ayi != 0.)
		{
			/* Computing MAX */
			tol = max(tol, atoli / ayi);
		}
	}
	L160:
	/* Computing MAX */
	tol = max(tol, numeric_limits<double>::epsilon() * 100.);
	tol = min(tol, .001);
	sum = VectorMaxNorm(dls001.n, &rwork[lf0], &rwork[dls001.lewt-1]);
	/* Computing 2nd power */
	sum = 1. / (tol * w0 * w0) + tol * (sum * sum);
	h0 = 1. / sqrt(sum);
	h0 = min(h0, tdist);
	h0 = ((tout - ODEs.GetTime()) < 0) ? -h0 : h0;
	L180:
	rh = fabs(h0) * dls001.hmxi;
	if (rh > 1.)
	{
		h0 /= rh;
	}
	dls001.h = h0;
	for (i = 0; i < dls001.n; ++i)
	{
		rwork[i + lf0] = h0 * rwork[i + lf0];
	}
	goto L270;
	L200:
	nslast = dls001.nst;
	switch (itask)
	{
		case 1:
			goto L210;
		case 2:
			goto L250;
		case 3:
			goto L220;
		case 4:
			goto L230;
		case 5:
			goto L240;
	}
	L210:
	if ((dls001.tn - tout) * dls001.h < 0.)
	{
		goto L250;
	}
	if (InterpolateVariables(tout) != 0)
	{
		goto L627;
	}
	ODEs.SetTime(tout);
	goto L420;
	L220:
	tp = dls001.tn - dls001.hu * (numeric_limits<double>::epsilon() * 100. + 1.);
	if ((tp - tout) * dls001.h > 0.)
	{
		goto L623;
	}
	if ((dls001.tn - tout) * dls001.h < 0.)
	{
		goto L250;
	}
	ODEs.SetTime(dls001.tn);
	goto L400;
	L230:
	tcrit = rwork[0];
	if ((dls001.tn - tcrit) * dls001.h > 0.)
	{
		goto L624;
	}
	if ((tcrit - tout) * dls001.h < 0.)
	{
		goto L625;
	}
	if ((dls001.tn - tout) * dls001.h < 0.)
	{
		goto L245;
	}
	if (InterpolateVariables(tout) != 0)
	{
		goto L627;
	}
	ODEs.SetTime(tout);
	goto L420;
	L240:
	tcrit = rwork[0];
	if ((dls001.tn - tcrit) * dls001.h > 0.)
	{
		goto L624;
	}
	L245:
	hmx = fabs(dls001.tn) + fabs(dls001.h);
	ihit = fabs(dls001.tn - tcrit) <= numeric_limits<double>::epsilon() * 100. * hmx;
	if (ihit)
	{
		ODEs.SetTime(tcrit);
	}
	if (ihit)
	{
		goto L400;
	}
	tnext = dls001.tn + dls001.h * (numeric_limits<double>::epsilon() * 4. + 1.);
	if ((tnext - tcrit) * dls001.h <= 0.)
	{
		goto L250;
	}
	dls001.h = (tcrit - dls001.tn) * (1. - numeric_limits<double>::epsilon() * 4.);
	if (istate == 2 && dls001.jstart >= 0)
	{
		dls001.jstart = -2;
	}
	L250:
	if (dls001.meth == dlsa01.mused)
	{
		goto L255;
	}
	if (insufr == 1)
	{
		goto L550;
	}
	if (insufi == 1)
	{
		goto L555;
	}
	L255:
	if (dls001.nst - nslast >= mxstep)
	{
		goto L500;
	}
	ODEs.GetErrorWeights(itol, rtol, atol, &rwork[20], &rwork[dls001.lewt-1]);
	for (i = 0; i < dls001.n; ++i)
	{
		if (rwork[i + dls001.lewt - 1] <= 0.)
		{
			goto L510;
		}
		rwork[i + dls001.lewt - 1] = 1. / rwork[i + dls001.lewt - 1];
	}
	L270:
	tolsf = numeric_limits<double>::epsilon() * VectorMaxNorm(dls001.n, &rwork[20], &rwork[dls001.lewt-1]);
	if (tolsf <= 1.)
	{
		goto L280;
	}
	tolsf *= 2.;
	if (dls001.nst == 0)
	{
		goto L626;
	}
	goto L520;
	L280:
	if (dls001.tn + dls001.h != dls001.tn)
	{
		goto L290;
	}
	++nhnil;
	if (nhnil > mxhnil)
	{
		goto L290;
	}
	Message = "DLSODA - Warning..Internal T (=R1) and H (=R2) are such that in the machine, T + H = T on the next step (H = step size). Solver will continue anyway.";
	em.Message(Message, 101, 0, 0, 0, 0, 2, dls001.tn, dls001.h);
	if (nhnil < mxhnil)
	{
		goto L290;
	}
	Message = "DLSODA - Above warning has been issued I1 times. It will not be issued again for this problem.";
	em.Message(Message, 102, 0, 1, mxhnil, 0, 0, 0., 0.);
	L290:
	Step();
	kgo = 1 - dls001.kflag;
	switch (kgo)
	{
		case 1:
			goto L300;
		case 2:
			goto L530;
		case 3:
			goto L540;
	}
	L300:
	init = 1;
	if (dls001.meth == dlsa01.mused)
	{
		goto L310;
	}
	tsw = dls001.tn;
	dls001.maxord = dlsa01.mxordn;
	if (dls001.meth == 2)
	{
		dls001.maxord = dlsa01.mxords;
	}
	if (dls001.meth == 2)
	{
		rwork[dls001.lwm-1] = sqrt(numeric_limits<double>::epsilon());
	}
	insufr = min(insufr, (long)1);
	insufi = min(insufi, (long)1);
	dls001.jstart = -1;
	if (ixpr == 0)
	{
		goto L310;
	}
	if (dls001.meth == 2)
	{
		Message = "DLSODA - A switch to the BDF (stiff) method has occurred";
		em.Message(Message, 105, 0, 0, 0, 0, 0, 0., 0.);
	}
	if (dls001.meth == 1)
	{
		Message = "DLSODA - A switch to the Adams (nonstiff) method has occurred";
		em.Message(Message, 106, 0, 0, 0, 0, 0, 0., 0.);
	}
	Message = "at T = R1, tentative step size H = R2, step NST = I1";
	em.Message(Message, 107, 0, 1, dls001.nst, 0, 2, dls001.tn, dls001.h);
	L310:
	switch (itask)
	{
		case 1:
			goto L320;
		case 2:
			goto L400;
		case 3:
			goto L330;
		case 4:
			goto L340;
		case 5:
			goto L350;
	}
	L320:
	if ((dls001.tn - tout) * dls001.h < 0.)
	{
		goto L250;
	}
	InterpolateVariables(tout);
	ODEs.SetTime(tout);
	goto L420;
	L330:
	if ((dls001.tn - tout) * dls001.h >= 0.)
	{
		goto L400;
	}
	goto L250;
	L340:
	if ((dls001.tn - tout) * dls001.h < 0.)
	{
		goto L345;
	}
	InterpolateVariables(tout);
	ODEs.SetTime(tout);
	goto L420;
	L345:
	hmx = fabs(dls001.tn) + fabs(dls001.h);
	ihit = fabs(dls001.tn - tcrit) <= numeric_limits<double>::epsilon() * 100. * hmx;
	if (ihit)
	{
		goto L400;
	}
	tnext = dls001.tn + dls001.h * (numeric_limits<double>::epsilon() * 4. + 1.);
	if ((tnext - tcrit) * dls001.h <= 0.)
	{
		goto L250;
	}
	dls001.h = (tcrit - dls001.tn) * (1. - numeric_limits<double>::epsilon() * 4.);
	if (dls001.jstart >= 0)
	{
		dls001.jstart = -2;
	}
	goto L250;
	L350:
	hmx = fabs(dls001.tn) + fabs(dls001.h);
	ihit = fabs(dls001.tn - tcrit) <= numeric_limits<double>::epsilon() * 100. * hmx;
	L400:
	for (i = 0; i < dls001.n; ++i)
	{
		ODEs.SetVariable(i, rwork[i + 20]);
	}
	ODEs.SetTime(dls001.tn);
	if (itask != 4 && itask != 5)
	{
		goto L420;
	}
	if (ihit)
	{
		ODEs.SetTime(tcrit);
	}
	L420:
	istate = 2;
	rwork[10] = dls001.hu;
	rwork[11] = dls001.h;
	rwork[12] = dls001.tn;
	rwork[14] = tsw;
	iwork[10] = dls001.nst;
	iwork[13] = dls001.nqu;
	iwork[14] = dls001.nq;
	iwork[18] = dlsa01.mused;
	iwork[19] = dls001.meth;
	return;
	L500:
	Message = "DLSODA - At current T (=R1), MXSTEP (=I1) steps taken on this call before reaching TOUT";
	em.Message(Message, 201, 0, 1, mxstep, 0, 1, dls001.tn, 0.);
	istate = -1;
	goto L580;
	L510:
	ewti = rwork[dls001.lewt + i - 1];
	Message = "DLSODA - At T (=R1), EWT(I1) has become R2 .le. 0.";
	em.Message(Message, 202, 0, 1, i, 0, 2, dls001.tn, ewti);
	istate = -6;
	goto L580;
	L520:
	Message = "DLSODA - At T (=R1), too much accuracy requested for precision of machine.. See TOLSF (=R2)";
	em.Message(Message, 203, 0, 0, 0, 0, 2, dls001.tn, tolsf);
	rwork[13] = tolsf;
	istate = -2;
	goto L580;
	L530:
	Message = "DLSODA - At T(=R1) and step size H(=R2), the error test failed repeatedly or with ABS(H) = HMIN";
	em.Message(Message, 204, 0, 0, 0, 0, 2, dls001.tn, dls001.h);
	istate = -4;
	goto L560;
	L540:
	Message = "DLSODA - At T (=R1) and step size H (=R2), the corrector convergence failed repeatedly or with ABS(H) = HMIN";
	em.Message(Message, 205, 0, 0, 0, 0, 2, dls001.tn, dls001.h);
	istate = -5;
	goto L560;
	L550:
	Message = "DLSODA - At current T(=R1), RWORK length too small to proceed. The integration was otherwise successful.";
	em.Message(Message, 206, 0, 0, 0, 0, 1, dls001.tn, 0.);
	istate = -7;
	goto L580;
	L555:
	Message = "DLSODA - At current T(=R1), IWORK length too small to proceed. The integration was otherwise successful.";
	em.Message(Message, 207, 0, 0, 0, 0, 1, dls001.tn, 0.);
	istate = -7;
	goto L580;
	L560:
	big = 0.;
	imxer = 1;
	for (i = 0; i < dls001.n; ++i)
	{
		size = fabs(rwork[i + dls001.lacor - 1] * rwork[i + dls001.lewt - 1]);
		if (big < size)
		{
			big = size;
			imxer = i+1;
		}
	}
	iwork[15] = imxer;
	L580:
	for (i = 0; i < dls001.n; ++i)
	{
		ODEs.SetVariable(i, rwork[i + 20]);
	}
	ODEs.SetTime(dls001.tn);
	rwork[10] = dls001.hu;
	rwork[11] = dls001.h;
	rwork[12] = dls001.tn;
	rwork[14] = tsw;
	iwork[10] = dls001.nst;
	iwork[13] = dls001.nqu;
	iwork[14] = dls001.nq;
	iwork[18] = dlsa01.mused;
	iwork[19] = dls001.meth;
	return;
	L601:
	Message = "DLSODA - ISTATE (=I1) illegal.";
	em.Message(Message, 1, 0, 1, istate, 0, 0, 0., 0.);
	if (istate < 0)
	{
		goto L800;
	}
	goto L700;
	L602:
	Message = "DLSODA - ITASK (=I1) illegal.";
	em.Message(Message, 2, 0, 1, itask, 0, 0, 0., 0.);
	goto L700;
	L603:
	Message = "DLSODA - ISTATE .gt. 1 but DLSODA not initialized.";
	em.Message(Message, 3, 0, 0, 0, 0, 0, 0., 0.);
	goto L700;
	//L604:
	//Message = "DLSODA - NEQ (=I1) .lt. 1";
	//em.Message(Message, 4, 0, 1, ODEs.nVariables(), 0, 0, 0., 0.);
	//goto L700;
	L605:
	Message = "DLSODA - ISTATE = 3 and NEQ increased (I1 to I2).";
	em.Message(Message, 5, 0, 2, dls001.n, ODEs.nVariables(), 0, 0., 0.);
	goto L700;
	L606:
	Message = "DLSODA - ITOL (=I1) illegal.";
	em.Message(Message, 6, 0, 1, itol, 0, 0, 0., 0.);
	goto L700;
	L607:
	Message = "DLSODA - IOPT (=I1) illegal.";
	em.Message(Message, 7, 0, 1, iopt, 0, 0, 0., 0.);
	goto L700;
	L608:
	Message = "DLSODA - JT (=I1) illegal.";
	em.Message(Message, 8, 0, 1, dlsa01.jtyp, 0, 0, 0., 0.);
	goto L700;
	L609:
	Message = "DLSODA - ML (=I1) illegal: .lt.0 or .ge.NEQ (=I2)";
	em.Message(Message, 9, 0, 2, ml, ODEs.nVariables(), 0, 0., 0.);
	goto L700;
	L610:
	Message = "DLSODA - MU (=I1) illegal: .lt.0 or .ge.NEQ (=I2)";
	em.Message(Message, 10, 0, 2, mu, ODEs.nVariables(), 0, 0., 0.);
	goto L700;
	L611:
	Message = "DLSODA - IXPR (=I1) illegal.";
	em.Message(Message, 11, 0, 1, ixpr, 0, 0, 0., 0.);
	goto L700;
	L612:
	Message = "DLSODA - MXSTEP (=I1) .lt. 0";
	em.Message(Message, 12, 0, 1, mxstep, 0, 0, 0., 0.);
	goto L700;
	L613:
	Message = "DLSODA - MXHNIL (=I1) .lt. 0";
	em.Message(Message, 13, 0, 1, mxhnil, 0, 0, 0., 0.);
	goto L700;
	L614:
	Message = "DLSODA - TOUT (=R1) behind T (=R2)";
	em.Message(Message, 14, 0, 0, 0, 0, 2, tout, ODEs.GetTime());
	Message = "Integration direction is given by H0 (=R1)";
	em.Message(Message, 14, 0, 0, 0, 0, 1, h0, 0.);
	goto L700;
	L615:
	Message = "DLSODA - HMAX (=R1) .lt. 0.0";
	em.Message(Message, 15, 0, 0, 0, 0, 1, hmax, 0.);
	goto L700;
	L616:
	Message = "DLSODA - HMIN (=R1) .lt. 0.0";
	em.Message(Message, 16, 0, 0, 0, 0, 1, dls001.hmin, 0.);
	goto L700;
	L617:
	Message = "DLSODA - RWORK length needed, LENRW (=I1), exceeds LRW (=I2)";
	em.Message(Message, 17, 0, 2, lenrw, lrw, 0, 0., 0.);
	goto L700;
	L618:
	Message = "DLSODA - IWORK length needed, LENIW (=I1), exceeds LIW (=I2)";
	em.Message(Message, 18, 0, 2, leniw, liw, 0, 0., 0.);
	goto L700;
	L619:
	Message = "DLSODA - RTOL(I1) is R1 .lt. 0.0";
	em.Message(Message, 19, 0, 1, i, 0, 1, rtoli, 0.);
	goto L700;
	L620:
	Message = "DLSODA - ATOL(I1) is R1 .lt. 0.0";
	em.Message(Message, 20, 0, 1, i, 0, 1, atoli, 0.);
	goto L700;
	L621:
	ewti = rwork[dls001.lewt + i - 1];
	Message = "DLSODA - EWT(I1) is R1 .le. 0.0";
	em.Message(Message, 21, 0, 1, i, 0, 1, ewti, 0.);
	goto L700;
	L622:
	Message = "DLSODA - TOUT(=R1) too close to T(=R2) to start integration.";
	em.Message(Message, 22, 0, 0, 0, 0, 2, tout, ODEs.GetTime());
	goto L700;
	L623:
	Message = "DLSODA - ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)";
	em.Message(Message, 23, 0, 1, itask, 0, 2, tout, tp);
	goto L700;
	L624:
	Message = "DLSODA - ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)";
	em.Message(Message, 24, 0, 0, 0, 0, 2, tcrit, dls001.tn);
	goto L700;
	L625:
	Message = "DLSODA - ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)";
	em.Message(Message, 25, 0, 0, 0, 0, 2, tcrit, tout);
	goto L700;
	L626:
	Message = "DLSODA - At start of problem, too much accuracy requested for precision of machine.. See TOLSF (=R1)";
	em.Message(Message, 26, 0, 0, 0, 0, 1, tolsf, 0.);
	rwork[13] = tolsf;
	goto L700;
	L627:
	Message = "DLSODA - Trouble in DINTDY. ITASK = I1, TOUT = R1";
	em.Message(Message, 27, 0, 1, itask, 0, 1, tout, 0.);
	goto L700;
	L628:
	Message = "DLSODA - MXORDN (=I1) .lt. 0";
	em.Message(Message, 28, 0, 1, dlsa01.mxordn, 0, 0, 0., 0.);
	goto L700;
	L629:
	Message = "DLSODA - MXORDS (=I1) .lt. 0";
	em.Message(Message, 29, 0, 1, dlsa01.mxords, 0, 0, 0., 0.);
	L700:
	istate = -3;
	return;
	L800:
	Message = "DLSODA - Run aborted.. apparent infinite loop.";
	em.Message(Message, 303, 2, 0, 0, 0, 0, 0., 0.);
	return;
}

int LivermoreSolver::Solve(const double &time)
{
	dlsoda_(time);
	return istate;
}

