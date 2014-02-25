/*!  @file mc.cpp */

#include "Parseur.hpp"
#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"

void save(Container const& param, Container const& result, Write& w);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int nthreads(1);
	P.get("nthreads",nthreads);
	Container input;
	if(P.check("t_max")){
		input.set("t_max",P.get<unsigned int>("t_max"));
	}

	CreateSystem CS(P);
	Container param(true);
	CS.get_input(input);
	CS.get_param(param);
	Write file("bla.dat");

	if( CS.use_complex() ){
#pragma omp parallel num_threads(nthreads)
		{
			MonteCarlo<std::complex<double> > sim(input);
			sim.run(1,1e5);
			Container result(true);
			sim.save(result);
			save(param,result,file);
		}
	} else {
#pragma omp parallel num_threads(nthreads)
		{
			MonteCarlo<double> sim(input);
			sim.run(1,1e5);
			Container result(true);
			sim.save(result);
			save(param,result,file);
		}
	}
}

void save(Container const& param, Container const& result, Write& w){
	w<<"%";
	for(unsigned int i(0);i<param.size();i++){
		w<<param.name(i)<<" ";
	}
	for(unsigned int i(0);i<result.size();i++){
		w<<result.name(i)<<" ";
	}
	w<<Write::endl;
	w<<param<<" "<<result<<Write::endl;
}
