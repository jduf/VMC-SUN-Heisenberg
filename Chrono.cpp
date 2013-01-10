#include"Chrono.hpp"

// définitions des constructeurs
Chrono::Chrono()
	:tics(0),t(0.0)
{ }

void Chrono::tic(){
	tics = std::clock();
}

void Chrono::tac(){
	tics = (std::clock()-t);
	t = tics/double(CLOCKS_PER_SEC);
}

double const& Chrono::Get_t() const{
	return t;
}

// surcharge d'opérateurs externes
std::ostream& operator<<(std::ostream& flux, Chrono const& t){
	flux << t.Get_t();
	return flux;
}


