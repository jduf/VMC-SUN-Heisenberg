/*!  @file analyse.cpp */

#include "Read.hpp"
#include "Parseur.hpp"

void analyse_sim(std::string sim_name);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	analyse_sim(P.get<std::string>("0"));
}

void analyse_sim(std::string sim_name){
	unsigned int nruns;
	unsigned int N;
	unsigned int m;
	unsigned int n;
	int bc;
	double param;
	double E;
	double DeltaE;
	unsigned int Nsteps;
	unsigned int status;
	Vector<double> corr;
	Vector<unsigned int> ref;

	Read r(sim_name);
	r>>nruns>>ref>>N>>m>>n>>bc>>param;
	for(unsigned int j(0);j<nruns;j++){
		r>>E>>DeltaE>>Nsteps>>status>>corr;
		std::cout<<corr<<std::endl;
	}
	r>>E>>DeltaE>>corr;
	std::cout<<corr<<std::endl;
}
