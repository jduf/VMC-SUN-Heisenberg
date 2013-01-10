#ifndef DEF_CHRONO
#define DEF_CHRONO

#include<iostream>
#include<ctime>

class Chrono{
	private:
		clock_t tics;
		double t;

	public:
// définitions des constructeurs
		Chrono();

// accesseur
		double const& Get_t() const;
// methode de classe
		void tic(); // lance le chrono
		void tac(); // stoppe le chrono
};

// surcharge d'opérateurs externes
std::ostream& operator<<(std::ostream& flux, Chrono const& t);
#endif
