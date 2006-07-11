//
// Copyright (c) 2003-2005 João Abecasis
//

#ifndef _agedo_PARAMETERS_LIMITS_H
#define _agedo_PARAMETERS_LIMITS_H

#include "Resultset.h"

#include <limits>
#include <string>
#include <vector>

class Limits
{
public:
    unsigned int nIndividuals;
    unsigned int nGenerations;
    const unsigned int nParameters;
    double *Min, *Max;
    bool *AbsoluteMin, *AbsoluteMax;
    std::vector<ResultSet> ResultSets;
    std::string Best_name;

    Limits(unsigned int nParams, unsigned int nResultSets = 0) : nParameters(nParams), Min(0), Max(0), AbsoluteMin(0), AbsoluteMax(0), ResultSets(nResultSets)
    {
        nIndividuals = nParameters * 20;
        nGenerations = nParameters * 100;

        if (nParameters != 0)
        {
            AbsoluteMin = new bool[nParams];
            AbsoluteMax = new bool[nParams];
            Min = new double[nParams];
            Max = new double[nParams];

            for (size_t i = 0; i < nParameters; ++i)
            {
                AbsoluteMin[i] = false;
                AbsoluteMax[i] = false;
                Min[i] = std::numeric_limits<double>::min();
                Max[i] = std::numeric_limits<double>::max();
            }
        }
    }

    ~Limits()
    {
        delete [] AbsoluteMin;
        delete [] AbsoluteMax;
        delete [] Min;
        delete [] Max;
    }
};

#endif
