#include"Swarm.hpp"
#include"Parseur.hpp"

int main(int argc,char* argv[]){
	Parseur P(argc,argv);
	Function f;
	Swarm s(100,164,P.get<double>("cg"),P.get<double>("cp"),f.f);
	do{ s.next_step(); }
	while (s.check());

	//s.path();
	s.result();
}

