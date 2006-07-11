//
// Copyright (c) 2003 Joao Abecasis
//

#ifndef _agedo_LINEAR_ALGEBRA_2_H
#define _agedo_LINEAR_ALGEBRA_2_H

#include <algorithm>

namespace LSODA { namespace LinearAlgebra {

	using namespace std;

void dgefa_(double *a, const long lda, const long n, long *ipvt, long *info)
{
    /* Local variables */
    long j, k, l, kp1, nm1;
    double t;

/* ***BEGIN PROLOGUE  DGEFA */
/* ***PURPOSE  Factor a matrix using Gaussian elimination. */
/* ***CATEGORY  D2A1 */
/* ***TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C) */
/* ***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK, */
/*             MATRIX FACTORIZATION */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     DGEFA factors a double precision matrix by Gaussian elimination. */

/*     DGEFA is usually called by DGECO, but it can be called */
/*     directly with a saving in time if  RCOND  is not needed. */
/*     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) . */

/*     On Entry */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                the matrix to be factored. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*     On Return */

/*        A       an upper triangular matrix and the multipliers */
/*                which were used to obtain it. */
/*                The factorization can be written  A = L*U  where */
/*                L  is a product of permutation and unit lower */
/*                triangular matrices and  U  is upper triangular. */

/*        IPVT    INTEGER(N) */
/*                an integer vector of pivot indices. */

/*        INFO    INTEGER */
/*                = 0  normal value. */
/*                = K  if  U(K,K) .EQ. 0.0 .  This is not an error */
/*                     condition for this subroutine, but it does */
/*                     indicate that DGESL or DGEDI will divide by zero */
/*                     if called.  Use  RCOND  in DGECO for a reliable */
/*                     indication of singularity. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DGEFA */


/*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING */

    /* Parameter adjustments */
    a -= 1 + lda;
    --ipvt;

    /* Function Body */
    *info = 0;
    nm1 = n - 1;
    if (nm1 < 1) {
	goto L70;
    }
    for (k = 1; k <= nm1; ++k) {
	kp1 = k + 1;

/*        FIND L = PIVOT INDEX */

	l = IndexOfMax(n - k + 1, &a[k + k * lda]) + k - 1;
	ipvt[k] = l;

/*        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED */

	if (a[l + k * lda] == 0.) {
	    goto L40;
	}

/*           INTERCHANGE IF NECESSARY */

	if (l == k) {
	    goto L10;
	}
	t = a[l + k * lda];
	a[l + k * lda] = a[k + k * lda];
	a[k + k * lda] = t;
L10:

/*           COMPUTE MULTIPLIERS */

	t = -1. / a[k + k * lda];
	AX(n - k, t, &a[k + 1 + k * lda]);

/*           ROW ELIMINATION WITH COLUMN INDEXING */

	for (j = kp1; j <= n; ++j) {
	    t = a[l + j * lda];
	    if (l == k) {
		goto L20;
	    }
	    a[l + j * lda] = a[k + j * lda];
	    a[k + j * lda] = t;
L20:
	    AXplusY(n - k, t, &a[k + 1 + k * lda], &a[k + 1 + j * lda]);
/* L30: */
	}
	goto L50;
L40:
	*info = k;
L50:
/* L60: */
	;
    }
L70:
    ipvt[n] = n;
    if (a[n + n * lda] == 0.) {
	*info = n;
    }
    return;
}

void dgesl_(double *a, const long lda, const long n, long *ipvt, double *b, const long job)
{
    /* Local variables */
    long k, l, kb, nm1;
    double t;

/* ***BEGIN PROLOGUE  DGESL */
/* ***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the */
/*            factors computed by DGECO or DGEFA. */
/* ***CATEGORY  D2A1 */
/* ***TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     DGESL solves the double precision system */
/*     A * X = B  or  TRANS(A) * X = B */
/*     using the factors computed by DGECO or DGEFA. */

/*     On Entry */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                the output from DGECO or DGEFA. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                the order of the matrix  A . */

/*        IPVT    INTEGER(N) */
/*                the pivot vector from DGECO or DGEFA. */

/*        B       DOUBLE PRECISION(N) */
/*                the right hand side vector. */

/*        JOB     INTEGER */
/*                = 0         to solve  A*X = B , */
/*                = nonzero   to solve  TRANS(A)*X = B  where */
/*                            TRANS(A)  is the transpose. */

/*     On Return */

/*        B       the solution vector  X . */

/*     Error Condition */

/*        A division by zero will occur if the input factor contains a */
/*        zero on the diagonal.  Technically this indicates singularity */
/*        but it is often caused by improper arguments or improper */
/*        setting of LDA .  It will not occur if the subroutines are */
/*        called correctly and if DGECO has set RCOND .GT. 0.0 */
/*        or DGEFA has set INFO .EQ. 0 . */

/*     To compute  INVERSE(A) * C  where  C  is a matrix */
/*     with  P  columns */
/*           CALL DGECO(A,LDA,N,IPVT,RCOND,Z) */
/*           IF (RCOND is too small) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL DGESL(A,LDA,N,IPVT,C(1,J),0) */
/*        10 CONTINUE */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  DAXPY, DDOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DGESL */

    /* Parameter adjustments */
    a -= 1 + lda;
    --ipvt;
    --b;

    /* Function Body */
    nm1 = n - 1;
    if (job != 0) {
	goto L50;
    }

/*        JOB = 0 , SOLVE  A * X = B */
/*        FIRST SOLVE  L*Y = B */

    if (nm1 < 1) {
	goto L30;
    }
    for (k = 1; k <= nm1; ++k) {
	l = ipvt[k];
	t = b[l];
	if (l == k) {
	    goto L10;
	}
	b[l] = b[k];
	b[k] = t;
L10:
	AXplusY(n - k, t, &a[k + 1 + k * lda], &b[k + 1]);
/* L20: */
    }
L30:

/*        NOW SOLVE  U*X = Y */

    for (kb = 1; kb <= n; ++kb) {
	k = n + 1 - kb;
	b[k] /= a[k + k * lda];
	t = -b[k];
	AXplusY(k - 1, t, &a[k * lda + 1], &b[1]);
/* L40: */
    }
    goto L100;
L50:

/*        JOB = NONZERO, SOLVE  TRANS(A) * X = B */
/*        FIRST SOLVE  TRANS(U)*Y = B */

    for (k = 1; k <= n; ++k) {
	t = DotProduct(k - 1, &a[k * lda + 1], &b[1]);
	b[k] = (b[k] - t) / a[k + k * lda];
/* L60: */
    }

/*        NOW SOLVE TRANS(L)*X = Y */

    if (nm1 < 1) {
	goto L90;
    }
    for (kb = 1; kb <= nm1; ++kb) {
	k = n - kb;
	b[k] += DotProduct(n - k, &a[k + 1 + k * lda], &b[k + 1]);
	l = ipvt[k];
	if (l == k) {
	    goto L70;
	}
	t = b[l];
	b[l] = b[k];
	b[k] = t;
L70:
/* L80: */
	;
    }
L90:
L100:
    return;
}

void dgbfa_(double *abd, const long lda, const long n, long *ml, long *mu, long *ipvt, long *info)
{
    /* Local variables */
    long i__, j, k, l, m, i0, j0, j1, lm, mm, ju, jz, kp1, nm1;
    double t;

/* ***BEGIN PROLOGUE  DGBFA */
/* ***PURPOSE  Factor a band matrix using Gaussian elimination. */
/* ***CATEGORY  D2A2 */
/* ***TYPE      DOUBLE PRECISION (SGBFA-S, DGBFA-D, CGBFA-C) */
/* ***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     DGBFA factors a double precision band matrix by elimination. */

/*     DGBFA is usually called by DGBCO, but it can be called */
/*     directly with a saving in time if  RCOND  is not needed. */

/*     On Entry */

/*        ABD     DOUBLE PRECISION(LDA, N) */
/*                contains the matrix in band storage.  The columns */
/*                of the matrix are stored in the columns of  ABD  and */
/*                the diagonals of the matrix are stored in rows */
/*                ML+1 through 2*ML+MU+1 of  ABD . */
/*                See the comments below for details. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  ABD . */
/*                LDA must be .GE. 2*ML + MU + 1 . */

/*        N       INTEGER */
/*                the order of the original matrix. */

/*        ML      INTEGER */
/*                number of diagonals below the main diagonal. */
/*                0 .LE. ML .LT.  N . */

/*        MU      INTEGER */
/*                number of diagonals above the main diagonal. */
/*                0 .LE. MU .LT.  N . */
/*                More efficient if  ML .LE. MU . */
/*     On Return */

/*        ABD     an upper triangular matrix in band storage and */
/*                the multipliers which were used to obtain it. */
/*                The factorization can be written  A = L*U  where */
/*                L  is a product of permutation and unit lower */
/*                triangular matrices and  U  is upper triangular. */

/*        IPVT    INTEGER(N) */
/*                an integer vector of pivot indices. */

/*        INFO    INTEGER */
/*                = 0  normal value. */
/*                = K  if  U(K,K) .EQ. 0.0 .  This is not an error */
/*                     condition for this subroutine, but it does */
/*                     indicate that DGBSL will divide by zero if */
/*                     called.  Use  RCOND  in DGBCO for a reliable */
/*                     indication of singularity. */

/*     Band Storage */

/*           If  A  is a band matrix, the following program segment */
/*           will set up the input. */

/*                   ML = (band width below the diagonal) */
/*                   MU = (band width above the diagonal) */
/*                   M = ML + MU + 1 */
/*                   DO 20 J = 1, N */
/*                      I1 = MAX(1, J-MU) */
/*                      I2 = MIN(N, J+ML) */
/*                      DO 10 I = I1, I2 */
/*                         K = I - J + M */
/*                         ABD(K,J) = A(I,J) */
/*                10    CONTINUE */
/*                20 CONTINUE */

/*           This uses rows  ML+1  through  2*ML+MU+1  of  ABD . */
/*           In addition, the first  ML  rows in  ABD  are used for */
/*           elements generated during the triangularization. */
/*           The total number of rows needed in  ABD  is  2*ML+MU+1 . */
/*           The  ML+MU by ML+MU  upper left triangle and the */
/*           ML by ML  lower right triangle are not referenced. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DGBFA */


/* ***FIRST EXECUTABLE STATEMENT  DGBFA */
    /* Parameter adjustments */
    abd -= 1 + lda;
    --ipvt;

    /* Function Body */
    m = *ml + *mu + 1;
    *info = 0;

/*     ZERO INITIAL FILL-IN COLUMNS */

    j0 = *mu + 2;
    j1 = min(n, m) - 1;
    if (j1 < j0) {
	goto L30;
    }
    for (jz = j0; jz <= j1; ++jz) {
	i0 = m + 1 - jz;
	for (i__ = i0; i__ <= *ml; ++i__) {
	    abd[i__ + jz * lda] = 0.;
/* L10: */
	}
/* L20: */
    }
L30:
    jz = j1;
    ju = 0;

/*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING */

    nm1 = n - 1;
    if (nm1 < 1) {
	goto L130;
    }
    for (k = 1; k <= nm1; ++k) {
	kp1 = k + 1;

/*        ZERO NEXT FILL-IN COLUMN */

	++jz;
	if (jz > n) {
	    goto L50;
	}
	if (*ml < 1) {
	    goto L50;
	}
	for (i__ = 1; i__ <= *ml; ++i__) {
	    abd[i__ + jz * lda] = 0.;
/* L40: */
	}
L50:

/*        FIND L = PIVOT INDEX */

/* Computing MIN */
	lm = min(*ml, n - k);
	l = IndexOfMax(lm + 1, &abd[m + k * lda]) + m - 1;
	ipvt[k] = l + k - m;

/*        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED */

	if (abd[l + k * lda] == 0.) {
	    goto L100;
	}

/*           INTERCHANGE IF NECESSARY */

	if (l == m) {
	    goto L60;
	}
	t = abd[l + k * lda];
	abd[l + k * lda] = abd[m + k * lda];
	abd[m + k * lda] = t;
L60:

/*           COMPUTE MULTIPLIERS */

	t = -1. / abd[m + k * lda];
	AX(lm, t, &abd[m + 1 + k * lda]);

/*           ROW ELIMINATION WITH COLUMN INDEXING */

/* Computing MIN */
/* Computing MAX */
	ju = min(max(ju, *mu + ipvt[k]), n);
	mm = m;
	if (ju < kp1) {
	    goto L90;
	}
	for (j = kp1; j <= ju; ++j) {
	    --l;
	    --mm;
	    t = abd[l + j * lda];
	    if (l == mm) {
		goto L70;
	    }
	    abd[l + j * lda] = abd[mm + j * lda];
	    abd[mm + j * lda] = t;
L70:
	    AXplusY(lm, t, &abd[m + 1 + k * lda], &abd[mm + 1 + j * lda]);
/* L80: */
	}
L90:
	goto L110;
L100:
	*info = k;
L110:
/* L120: */
	;
    }
L130:
    ipvt[n] = n;
    if (abd[m + n * lda] == 0.) {
	*info = n;
    }
    return;
}

void dgbsl_(double *abd, const long lda, const long n, long *ml, long *mu, long *ipvt, double *b, const long job)
{
    /* Local variables */
    long k, l, m, kb, la, lb, lm, nm1;
    double t;

/* ***BEGIN PROLOGUE  DGBSL */
/* ***PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using */
/*            the factors computed by DGBCO or DGBFA. */
/* ***CATEGORY  D2A2 */
/* ***TYPE      DOUBLE PRECISION (SGBSL-S, DGBSL-D, CGBSL-C) */
/* ***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE */
/* ***AUTHOR  Moler, C. B., (U. of New Mexico) */
/* ***DESCRIPTION */

/*     DGBSL solves the double precision band system */
/*     A * X = B  or  TRANS(A) * X = B */
/*     using the factors computed by DGBCO or DGBFA. */

/*     On Entry */

/*        ABD     DOUBLE PRECISION(LDA, N) */
/*                the output from DGBCO or DGBFA. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  ABD . */

/*        N       INTEGER */
/*                the order of the original matrix. */

/*        ML      INTEGER */
/*                number of diagonals below the main diagonal. */

/*        MU      INTEGER */
/*                number of diagonals above the main diagonal. */

/*        IPVT    INTEGER(N) */
/*                the pivot vector from DGBCO or DGBFA. */

/*        B       DOUBLE PRECISION(N) */
/*                the right hand side vector. */

/*        JOB     INTEGER */
/*                = 0         to solve  A*X = B , */
/*                = nonzero   to solve  TRANS(A)*X = B , where */
/*                            TRANS(A)  is the transpose. */

/*     On Return */

/*        B       the solution vector  X . */

/*     Error Condition */

/*        A division by zero will occur if the input factor contains a */
/*        zero on the diagonal.  Technically this indicates singularity */
/*        but it is often caused by improper arguments or improper */
/*        setting of LDA .  It will not occur if the subroutines are */
/*        called correctly and if DGBCO has set RCOND .GT. 0.0 */
/*        or DGBFA has set INFO .EQ. 0 . */

/*     To compute  INVERSE(A) * C  where  C  is a matrix */
/*     with  P  columns */
/*           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z) */
/*           IF (RCOND is too small) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0) */
/*        10 CONTINUE */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  DAXPY, DDOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DGBSL */

/* ***FIRST EXECUTABLE STATEMENT  DGBSL */
    /* Parameter adjustments */
    abd -= 1 + lda;
    --ipvt;
    --b;

    /* Function Body */
    m = *mu + *ml + 1;
    nm1 = n - 1;
    if (job != 0) {
	goto L50;
    }

/*        JOB = 0 , SOLVE  A * X = B */
/*        FIRST SOLVE L*Y = B */

    if (*ml == 0) {
	goto L30;
    }
    if (nm1 < 1) {
	goto L30;
    }
    for (k = 1; k <= nm1; ++k) {
/* Computing MIN */
	lm = min(*ml, n - k);
	l = ipvt[k];
	t = b[l];
	if (l == k) {
	    goto L10;
	}
	b[l] = b[k];
	b[k] = t;
L10:
	AXplusY(lm, t, &abd[m + 1 + k * lda], &b[k + 1]);
/* L20: */
    }
L30:

/*        NOW SOLVE  U*X = Y */

    for (kb = 1; kb <= n; ++kb) {
	k = n + 1 - kb;
	b[k] /= abd[m + k * lda];
	lm = min(k,m) - 1;
	la = m - lm;
	lb = k - lm;
	t = -b[k];
	AXplusY(lm, t, &abd[la + k * lda], &b[lb]);
/* L40: */
    }
    goto L100;
L50:

/*        JOB = NONZERO, SOLVE  TRANS(A) * X = B */
/*        FIRST SOLVE  TRANS(U)*Y = B */

    for (k = 1; k <= n; ++k) {
	lm = min(k,m) - 1;
	la = m - lm;
	lb = k - lm;
	t = DotProduct(lm, &abd[la + k * lda], &b[lb]);
	b[k] = (b[k] - t) / abd[m + k * lda];
/* L60: */
    }

/*        NOW SOLVE TRANS(L)*X = Y */

    if (*ml == 0) {
	goto L90;
    }
    if (nm1 < 1) {
	goto L90;
    }
    for (kb = 1; kb <= nm1; ++kb) {
	k = n - kb;
/* Computing MIN */
	lm = min(*ml, n - k);
	b[k] += DotProduct(lm, &abd[m + 1 + k * lda], &b[k + 1]);
	l = ipvt[k];
	if (l == k) {
	    goto L70;
	}
	t = b[l];
	b[l] = b[k];
	b[k] = t;
L70:
/* L80: */
	;
    }
L90:
L100:
    return;
}

} /* namespace LSODA */ } /* namespace LinearAlgebra */

#endif
