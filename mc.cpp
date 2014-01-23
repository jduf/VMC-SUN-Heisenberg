/*!  @file mc.cpp */

#include "Parseur.hpp"
#include "MonteCarlo.hpp"
#include "ExtractSystem.hpp"

void save(Container const& param, Container const& output, Write& w);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::string filename(P.get<std::string>("sim"));
	unsigned int nthreads(1);
	if(!P.status()){
		Container input;
		Container param(true);
		ExtractSystem system(filename);
		input.set("filename",filename);
		Write results(filename+".dat");

		P.get("nthreads",nthreads);
		if(P.check("t_max")){
			input.set("t_max",P.get<unsigned int>("t_max"));
		}
		system.extract(input,param);

		if( system.use_complex() ){
#pragma omp parallel num_threads(nthreads)
			{
				MonteCarlo<std::complex<double> > sim(input);
				sim.run(1,1e5);
				save(param,sim.save(),results);
			}
		} else {
#pragma omp parallel num_threads(nthreads)
			{
				MonteCarlo<double> sim(input);
				sim.run(1,1e5);
				save(param,sim.save(),results);
			}
		}
	} else {
		fprintf(stderr,"main : ./mc -sim filename.jdbin -nthreads n\n");
	}
}

void save(Container const& param, Container const& output, Write& w){
	w<<"%";
	for(unsigned int i(0);i<param.size();i++){
		w<<param.name(i)<<" ";
	}
	for(unsigned int i(0);i<output.size();i++){
		w<<output.name(i)<<" ";
	}
	w<<Write::endl;
	w<<param<<" "<<output<<Write::endl;
}
