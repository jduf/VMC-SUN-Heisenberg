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
		input.set("filename",filename);

		std::string wf;
		unsigned int m(0),N(0),nthreads(1);
		bool fermionic(true);

		P.get("nthreads",nthreads);

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

		Write results(filename+".dat");


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
			if( wf == "jastrow" || wf == "trianglejastrow" ){ 
				Vector<double> nu;
				file.extract<Vector<double> >(nu);
				file.extract<Matrix<unsigned int> >("nn",input);
				file.extract<Vector<unsigned int> >("sl",input);
				file.extract<Matrix<std::complex<double> > >("omega",input);
				fermionic=false;
				input.set("nu",nu);
				param.set("nu",nu);
			}
		}
		input.set("fermionic",fermionic);
		file.extract<Matrix<unsigned int> >("sts",input);

		if( wf != "csl" && wf != "phi" && wf != "trianglephi" && wf != "jastrow" && wf != "trianglejastrow"){
			if(fermionic){
				file.extract<Matrix<double> >("EVec",input);
			}
#pragma omp parallel num_threads(nthreads)
			{
				MonteCarlo<double> sim(input);
				sim.run(1,1e5);
				save(param,sim.save(),results);
			}
		} else {
			std::cerr<<"complex"<<std::endl;
			if(fermionic){
				file.extract<Matrix<std::complex<double> > >("EVec",input);
			}
#pragma omp parallel num_threads(nthreads)
			{
				MonteCarlo<std::complex<double> > sim(input);
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
