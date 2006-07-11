//
// Copyright (c) 2004 Joao Abecasis / António Ferreira
//

#ifndef _agedo_RESULTSET_H
#define _agedo_RESULTSET_H

#include <deque>
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

	unsigned int nVariables;
	double range;
	bool discard;
	deque<ResultValue> Results;

	ResultSet()
	{
		nVariables = 0;
		discard = false;
		range = 0.0;
	}

	unsigned int ReadFromFile(const char * filename, unsigned int StartAtLine = 0, unsigned int EndAtLine = -1)
	// Returns number of lines actually read from file, including lines skipped at the beggining of file
	{
		ifstream data(filename);
		string strLine;

		unsigned int i, j, pos;
		const char * pcCh;
		char * pCh;
		ResultValue rv;
		double result;

		for (i = 0; i < StartAtLine; ++i)
			getline(data, strLine);

		while (i <= EndAtLine && !data.eof()){
			getline(data, strLine);
			pos = 0;
			pos = strLine.find_first_not_of(" \t",pos);
	    if (pos == string::npos) continue;

      rv.NotANumber = !getValidNumber(strLine, result, pos);
      rv.Value = result;
      if (rv.NotANumber) continue;  // skip lines that not start with a valid number

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

	void computeRange()
	// Returns number of lines actually read from file, including lines skipped at the beggining of file
	{
		double maxvalue = 0.0;
		double minvalue = 0.0;
		unsigned int j, k;
		double value;

		// Look in initial values
		for (j = 0; j < nVariables; ++j)
		{
			value = Results[j + 1].Value; // time is skipped
			if (value < minvalue) minvalue = value;
			if (value > maxvalue) maxvalue = value;
		}

		for (j = nVariables + 1; j < Results.size(); j += nVariables + 1)
		{
			for (k = 0; k < nVariables; ++k)
			{
				if (!Results[j + k + 1].NotANumber)
				{
					value = Results[j + k + 1].Value; // time is skipped
					if (value < minvalue) minvalue = value;
					if (value > maxvalue) maxvalue = value;
				}
			}
		}
		range = abs(maxvalue-minvalue);
	}

	private:
	bool getValidNumber(const string & input, double & output, unsigned & pos) {
	     // skip whitespace
	    if (pos >= input.size()) return false;
	     pos = input.find_first_not_of(" \t",pos);
	     if (pos == string::npos) return false;

	     string prov = input.substr(pos, input.size());
	     const char * currinput = prov.c_str();
	     char * endchar;

	     double outprov = strtod(currinput,&endchar);
	     if (endchar != currinput){
	        output = outprov;
	        pos += (endchar-currinput);
	        return true;
	     }
	     else {
	       output = 0.0;
	       return false;
	     }
	}

	// TODO: define public interface
};

#endif
