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
	tics = (std::clock()-t);
	t = tics/double(CLOCKS_PER_SEC);
}

bool Chrono::time_limit_reached(double sec){
	clock_t limit(sec*CLOCKS_PER_SEC);
	tics = (std::clock()-t);
	if(tics >= limit){return true;}
	else{return false;}
}
/*}*/

std::ostream& operator<<(std::ostream& flux, Chrono const& t){
	flux << t.get_t();
	return flux;
}

