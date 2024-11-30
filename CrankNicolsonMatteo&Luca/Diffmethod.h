#ifndef DIFFMETHOD_H
#define DIFFMETHOD_H
#include "Matrix.h"
#include "Option.h"
#include "Utils.h"
#include <iostream>
#include <cmath>
#include <vector>

namespace m2
{
	class American : public Option
	{
	private:
		Matrix values_;
		unsigned int N_ = timeDiscr_;
		unsigned int M_ = spotDiscr_;

	public:
		American();
		
		void priceCall();
		void pricePut();
	};

	class European : public Option
	{
	private:
		Matrix values_;

	public:
		European();

		void priceCall();
		void pricePut();
	};
}

#endif // !DIFFMETHOD_H

