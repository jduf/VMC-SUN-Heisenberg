/*!  @file mc.cpp */

#include "MonteCarlo.hpp"

std::string init(CreateSystem const& cs);
template<typename Type>
void run(CreateSystem const& cs, std::string const& path, unsigned int const& nruns, unsigned int const& tmax, unsigned int const& type);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int type(P.get<unsigned int>("type"));
	unsigned int nruns(P.get<unsigned int>("nruns"));
	unsigned int tmax(P.get<unsigned int>("tmax"));
	CreateSystem cs(P);
	if(!P.status()){
		std::string path(init(cs));
		if(P.is_vector("param")){
			Vector<double> param(P.get<Vector<double> >("param"));
			for(unsigned int i(0);i<param.size();i++){
				cs.create(param(i),type);
				if(!cs.is_degenerate()){
					if(cs.use_complex()){ run<std::complex<double> >(cs,path,nruns,tmax,type); } 
					else { run<double>(cs,path,nruns,tmax,type); }
				}
			}
		} else {
			double param(P.get<double>("param"));
			cs.create(param,type);
			if(!cs.is_degenerate()){
				if(cs.use_complex()){ run<std::complex<double> >(cs,path,nruns,tmax,type); }
				else { run<double>(cs,path,nruns,tmax,type); }
			}
		}
	}
}

std::string init(CreateSystem const& cs){
	std::string path("N"+tostring(cs.get_system().get_N())+"/m"+tostring(cs.get_system().get_m())+"/n"+tostring(cs.get_system().get_n())+"/");
	switch(cs.get_system().get_bc()){
		case -1:{path += "A/";}break;
		case 0: {path += "O/";}break;
		case 1: {path += "P/";}break;
		default:{std::cerr<<"Unknown boundary condition"<<std::endl;}
	}
	Linux command;
	command("mkdir -p " + path);
	return path;
}

template<typename Type>
void run(CreateSystem const& cs, std::string const& path, unsigned int const& nruns, unsigned int const& tmax, unsigned int const& type){
	IOFiles results(path+cs.get_filename()+".jdbin",true);
	results("type of simulation",type);
	RST rst;
	rst.title("Input","-");
	results.add_to_header(rst.get());
	rst.set();
	cs.save(results);
	results("number of simulations runned",nruns);
	rst.title("Results","-");
	results.add_to_header(rst.get());

	Data<double> E;
	DataSet<double> corr;
	DataSet<double> long_range_corr;
	E.set_conv(true);
	std::cout<<"will need to correct that"<<std::endl;
	corr.set(cs.get_system().get_n());
	if(type == 2){
		long_range_corr.set(cs.get_system().get_n()/3);
	}

#pragma omp parallel for 
	for(unsigned int i=0;i<nruns;i++){
		MCSystem<Type>* S(NULL);
		if(cs.is_bosonic()){ S = new SystemBosonic<Type>(cs,type); } 
		else { S = new SystemFermionic<Type>(cs,type); }
		MonteCarlo<Type> sim(S,tmax);
		sim.run();
#pragma omp critical
		{
			E.add_sample((sim.get_system())->get_energy());
			corr.add_sample((sim.get_system())->get_corr());
			long_range_corr.add_sample((sim.get_system())->get_long_range_corr());
			(sim.get_system())->save(results);
		}
		delete S;
	}

	E.complete_analysis();
	corr.complete_analysis();
		long_range_corr.complete_analysis();

	rst.set();
	rst.title("Mean results","-");
	results.add_to_header(rst.get());
	results("energy per site",E);
	results("correlation on links",corr);
	results("long range correlation",long_range_corr);

}
