//
// Copyright (c) 2004-2005 João Abecasis / António Ferreira
//

#ifndef _agedo_RESULTSET_H
#define _agedo_RESULTSET_H

#include <vector>
#include <fstream>
#include <string>

class ResultSet
{
public:
    struct ResultValue
    {
        bool NotANumber;
        double Value;
        ResultValue() : NotANumber(true) { }
    };

    unsigned int nVariables; // includes time

    std::vector<double> min;
    std::vector<double> max;
    std::vector<double> sum;
    std::vector<double> count;

    std::vector<double> range;
    std::vector<double> average;

    std::string title;

    bool discard;
    std::vector<ResultValue> Results;

    ResultSet()
    {
        nVariables = 0;
        discard = false;
    }

    unsigned int ReadFromFile(
        const char * filename,
        unsigned int StartAtLine = 0,
        unsigned int EndAtLine = (unsigned int)-1)
    // Returns number of lines actually read from file, including lines skipped at the beggining of file
    {
        if (title.empty())
            title = filename;

        std::ifstream data(filename);
        std::string strLine;

        unsigned int i, j, pos;
        ResultValue rv;
        double result;

        for (i = 0; i < StartAtLine; ++i)
            getline(data, strLine);

        while (i <= EndAtLine && !data.eof()){
            getline(data, strLine);
            pos = 0;
            pos = strLine.find_first_not_of(" \t",pos);
            if (pos == std::string::npos) continue;


            // skip lines that do not start with a valid number
            rv.NotANumber = !getValidNumber(strLine, result, pos);
            rv.Value = result;
            if (rv.NotANumber) continue;

            Results.push_back(rv);
            pos = strLine.find_first_of(" \t",pos);


            for (j = 1; j < nVariables; ++j){
                rv.NotANumber = !getValidNumber(strLine, result, pos);
                rv.Value = result;
                Results.push_back(rv);
                pos = strLine.find_first_of(" \t",pos);
            }
            ++i;
        }

        data.close();

        return --i;
    }

    void computeStatistics()
    {
        min.resize(nVariables);
        max.resize(nVariables);
        sum.resize(nVariables);
        count.resize(nVariables);

        range.resize(nVariables);
        average.resize(nVariables);

        for (unsigned int i = 0; i < nVariables; ++i)
        {
            // Assume initial values are always Numbers!
            min[i] = Results[i].Value;
            max[i] = Results[i].Value;
            sum[i] = Results[i].Value;
            count[i] = 1;
        }

        for (unsigned int i = nVariables; i < Results.size(); i += nVariables)
        {
            for (unsigned j = 0; j < nVariables; ++j)
            {
                if (!Results[i + j].NotANumber)
                {
                    double const & cur_value = Results[i + j].Value;

                    if ( cur_value < min[j])
                        min[j] = cur_value;
                    else if (max[j] < cur_value)
                        max[j] = cur_value;

                    sum[j] += cur_value;
                    ++count[j];
                }
            }
        }

        for (unsigned int i = 0; i < nVariables; ++i)
        {
            range[i] = max[i] - min[i];
            average[i] = sum[i] / count[i];
        }
    }

private:
    bool getValidNumber(const std::string & input, double & output, unsigned & pos)
    {
         // skip whitespace
        if (pos >= input.size()) return false;
        pos = input.find_first_not_of(" \t",pos);
        if (pos == std::string::npos) return false;

        std::string prov = input.substr(pos, input.size());
        const char * currinput = prov.c_str();
        char * endchar;

        double outprov = strtod(currinput,&endchar);
        if (endchar != currinput)
        {
            output = outprov;
            pos += (endchar-currinput);
            return true;
        }
        else
        {
            output = 0.0;
            return false;
        }
    }

    // TODO: define public interface
};

#endif
