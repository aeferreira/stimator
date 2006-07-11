// $Header: //194.117.40.12/c\044/Documents\040and\040Settings/CVSNT/Repositories/Bioquimica\040Teorica/AGEDO/Source/ga/GAEvalData.h,v 1.1 2003/09/05 16:46:14 Administrator Exp $
/* ----------------------------------------------------------------------------
  eval.h
  mbwall 3dec95
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved

 DESCRIPTION:
  This is the basic interface for the object that contains evaluation data.  It
can be used with genomes and/or populations in combination with their 
respective evaluation methods.
---------------------------------------------------------------------------- */
#ifndef _ga_eval_h_
#define _ga_eval_h_

class GAEvalData {
public:
  GAEvalData() {}
  GAEvalData(const GAEvalData&) {}
  virtual ~GAEvalData() {}
  GAEvalData& operator=(const GAEvalData& orig) 
    { if(&orig != this) copy(orig); return *this; }
  virtual GAEvalData* clone() const =0;
  virtual void copy(const GAEvalData&) =0;
};

#endif
