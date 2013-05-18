#include "Chrono.hpp"

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
	t = double(std::clock()-tics)/double(CLOCKS_PER_SEC);
}

bool Chrono::time_limit_reached(double sec) const{
	clock_t limit(sec*CLOCKS_PER_SEC);
	if(std::clock()-tics >= limit){return true;}
	else{return false;}
}

double const& Chrono::get_t() { 
	t = double(std::clock()-tics)/double(CLOCKS_PER_SEC);
	return t;
}
/*}*/

std::ostream& operator<<(std::ostream& flux, Chrono& t){
	flux << t.get_t();
	return flux;
}

