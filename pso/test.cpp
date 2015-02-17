#include"FuncPSO.hpp"
#include"Parseur.hpp"
#include"Time.hpp"

int main(int argc,char* argv[]){
	Parseur P(argc,argv);
	Time t;
	std::cout<<t.elapsed()<<std::endl;

	FuncPSO s(P.get<unsigned int>("Nbees"),P.get<unsigned int>("Nfreedom"),P.get<double>("cg"),P.get<double>("cp"),P.get<unsigned int>("maxiter"));

	s.PSO_init();
	std::cout<<"initialization done in "<<t.elapsed()<<std::endl;
	t.set();
	s.PSO_run();
	std::cout<<"run done in "<<t.elapsed()<<std::endl;
	s.result();
}
