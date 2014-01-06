#include"FuncPSO.hpp"
#include"Parseur.hpp"

int main(int argc,char* argv[]){
	Parseur P(argc,argv);
	for(unsigned int i(0);i<100;i++){
		FuncPSO s(P.get<unsigned int>("Nbees"),P.get<unsigned int>("Nfreedom"),P.get<double>("cg"),P.get<double>("cp"));
		//do{ s.next_step(); }
		//while (!s.converged());

		s.PSO_init(40);
		s.PSO_run();
		s.result();
		s.PSO_init(800);
		s.PSO_run(false);
		s.result();
	}
}

