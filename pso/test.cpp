#include"Swarm.hpp"
#include"Parseur.hpp"
#include"Functions.hpp"

int main(int argc,char* argv[]){
	Parseur P(argc,argv);
	Function F;
	Swarm<Function> s(P.get<unsigned int>("Nparticle"),P.get<unsigned int>("Nparam"),P.get<double>("cg"),P.get<double>("cp"),F,&Function::f);
	//do{ s.next_step(); }
	//while (!s.converged());

	s.run();
	//s.path();
	//s.result();
}

