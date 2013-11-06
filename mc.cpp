/*!  @file mc.cpp */

#include "Parseur.hpp"
#include "MonteCarlo.hpp"

#include<omp.h>

void save(Container const& param, Container const& output, Write& w);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::string filename(P.get<std::string>("sim"));
	if(!P.status()){
		Container input;
		Container param(true);
		std::string wf;
		unsigned int m(0),N(0),nthreads(1);

		P.get("nthreads",nthreads);

		Write results(filename+".dat");

		FileParser file(filename);
		file.extract<std::string>(wf);
		file.extract<unsigned int>(N);
		file.extract<unsigned int>(m);

		input.set("N",N);
		input.set("m",m);
		input.set("n",N*m);
		if(P.check("t_max")){
			input.set("t_max",P.get<unsigned int>("t_max"));
		}

		param.set("N",N);
		param.set("m",m);
		param.set("n",N*m);

		file.extract<double>("bc",param);
		if( wf != "chain" ){
			file.extract<unsigned int>("Lx",param);
			file.extract<unsigned int>("Ly",param);
			if(wf == "mu" || wf == "trianglemu"){
				file.extract<double>("mu",param);
			}
			if(wf == "phi" || wf == "trianglephi"){
				file.extract<double>("phi",param);
			}
		}
		file.extract<Matrix<unsigned int> >("sts",input);

		if( wf != "csl" && wf != "phi" && wf != "trianglephi"){
			MonteCarlo<double> sim(filename,nthreads);
			file.extract<Matrix<double> >("EVec",input);
#pragma omp parallel num_threads(nthreads)
			{
				sim.init(input,omp_get_thread_num());
				sim.run(omp_get_thread_num());
			}
			for(unsigned int thread(0); thread<nthreads; thread++){
				save(param,sim.save(thread),results);
			}
		} else {
			MonteCarlo<std::complex<double> > sim(filename,nthreads);
			file.extract<Matrix<std::complex<double> > >("EVec",input);
#pragma omp parallel num_threads(nthreads)
			{
				sim.init(input,omp_get_thread_num());
				sim.run(omp_get_thread_num());
			}
			for(unsigned int thread(0); thread<nthreads; thread++){
				save(param,sim.save(thread),results);
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
