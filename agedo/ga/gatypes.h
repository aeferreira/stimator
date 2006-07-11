// $Header: //194.117.40.12/c\044/Documents\040and\040Settings/CVSNT/Repositories/Bioquimica\040Teorica/AGEDO/Source/ga/gatypes.h,v 1.1 2003/09/05 16:46:17 Administrator Exp $
/* ----------------------------------------------------------------------------
  gatypes.h
  mbwall 29apr95
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved
---------------------------------------------------------------------------- */
#ifndef _ga_types_h_
#define _ga_types_h_


typedef float GAProbability;
typedef float GAProb;

// This enforces return values and arguments to functions that want strict
// boolean arguments.
typedef enum _GABoolean {
  gaFalse=0,
  gaTrue
} GABoolean, GABool;

typedef enum _GAStatus {
  gaSuccess=0,
  gaFailure=-1
} GAStatus;

// This defines what type to use for the bitstream data storage.  Use the
// smallest type available on your platform.  
typedef unsigned char GABit;


#define GAMax(a,b) ((a) > (b) ? (a) : (b))
#define GAMin(a,b) ((a) < (b) ? (a) : (b))


#endif
