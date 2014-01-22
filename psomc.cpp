/*!  @file psomc.cpp */

#include "PSOMonteCarlo.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	
	unsigned int Nfreedom(0);
	std::string wf(P.get<std::string>("wf"));
	unsigned int N(P.get<unsigned int>("N"));

	//if( wf == "jastrow" ){ Nfreedom=1; }
	//if( wf == "trianglejastrow" ){ Nfreedom=4; }
	if(N==2){ Nfreedom=1; }
	if(N==3){ Nfreedom=4; }

	PSOMonteCarlo s(P,Nfreedom);

	for(unsigned int i(0);i<Nfreedom;i++){
		s.PSO_set_limit(i,-1.5,1.5);
	}
	s.PSO_init();
	s.PSO_run(false); /*false because each run can vary in time*/
	s.save();
}
