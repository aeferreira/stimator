//
// Copyright (c) 2003 Joao Abecasis
//

#ifndef _agedo_ERROR_MESSENGER_H
#define _agedo_ERROR_MESSENGER_H

#include <iostream>
#include <string>

class ErrorMessenger
{

public:
	std::ostream* OutputStream;
	bool OutputMessages;
	unsigned long ExitLevel;

	ErrorMessenger() : OutputStream(&std::cerr), OutputMessages(true), ExitLevel(2)
	{
	}

	void Message(std::string msg, unsigned long nerr, unsigned long level, unsigned long ni, long i1, long i2, unsigned long nr, long double r1, long double r2)
	{
		if (OutputMessages)
		{
			*OutputStream << msg << std::endl;
			switch (ni)
			{
			case 1:
				*OutputStream << "In above message, I1 = " << i1 << std::endl;
				break;
			case 2:
				*OutputStream << "In above message, I1 = " << i1 << "; I2 = " << i2 << std::endl;
				break;
			}

			switch (nr)
			{
			case 1:
				*OutputStream << "In above message, R1 = " << r1 << std::endl;
				break;
			case 2:
				*OutputStream << "In above message, R1 = " << r1 << "; R2 = " << r2 << std::endl;
				break;
			}
		}

		if (level >= ExitLevel)
			exit(EXIT_FAILURE);

		return;
	}
};

#endif
