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

// accesseur (retourne le temps en seconde)
		double const& get_time() const;
// methode de classe
		void tic(); // lance le chrono
		void tac(); // stoppe le chrono
};

// surcharge d'opérateurs externes (retourne le temps en seconde)
std::ostream& operator<<(std::ostream& flux, Chrono const& t);

/*constructor*/
/*{*/
Chrono::Chrono()
	:tics(0),t(0.0)
{ }
/*}*/
/*methods*/
/*{*/
void Chrono::tic(){
	tics = std::clock();
}

void Chrono::tac(){
	tics = (std::clock()-t);
	t = tics/double(CLOCKS_PER_SEC);
}

double const& Chrono::get_time() const{
	return t;
}
/*}*/
/*operators*/
/*{*/
std::ostream& operator<<(std::ostream& flux, Chrono const& t){
	flux << t.get_time();
	return flux;
}
/*}*/
#endif
