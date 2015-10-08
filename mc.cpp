/*!  @file mc.cpp */

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"
#include <omp.h>

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int tmax(P.get<unsigned int>("tmax"));
	unsigned int i(0);
	unsigned int nruns(P.find("nruns",i,false)?P.get<unsigned int>(i):omp_get_max_threads());
	if(!P.find("M",i,false)){
		std::vector<unsigned int> M(P.get<unsigned int>("N"),P.get<unsigned int>("n")*P.get<unsigned int>("m")/P.get<unsigned int>("N"));
		P.set("M",M);
	}
	System sys(P);
	CreateSystem cs(&sys);
	cs.init(NULL,&P);
	if(!P.locked()){
		if(cs.get_status()==2){
			cs.create();
			if(cs.get_status()==1){
				cs.set_observables(P.find("ncorr",i,false)?P.get<unsigned int>(i):1);
				sys.set(cs.get_GS());
#pragma omp parallel for
				for(unsigned int j=0;j<nruns;j++){
					MCSystem* mcsys(NULL);
					if( cs.use_complex()){
						if(cs.is_bosonic()){
							mcsys = new SystemBosonic<std::complex<double> >(*dynamic_cast<const Bosonic<std::complex<double> >*>(cs.get_GS()));
						} else {
							mcsys = new SystemFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs.get_GS()));
						}
					} else {
						if(cs.is_bosonic()){
							mcsys = new SystemBosonic<double>(*dynamic_cast<const Bosonic<double>*>(cs.get_GS()));
						} else {
							mcsys = new SystemFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs.get_GS()));
						}
					}

					MonteCarlo sim(mcsys,tmax);
					sim.thermalize(1e6);
					sim.run();

#pragma omp critical
					{ sys.merge(mcsys); }

					delete mcsys;
				}

				sys.complete_analysis(1e-5);
				sys.delete_binning();

				std::cout<<sys.get_energy()<<std::endl;

				Linux command;
				command.mkpath(cs.get_path().c_str());
				IOFiles out(cs.get_path() + cs.get_filename()+".jdbin",true);
				cs.save_param(out);
				cs.save_input(out);
				sys.save_output(out);

				if(P.find("d",i,false)){
					CreateSystem tmp(&sys);
					tmp.init(NULL,&P);
					tmp.lattice("/tmp/",cs.get_filename());
					RSTFile rst("/tmp/",cs.get_filename());
					rst.figure("/tmp/"+cs.get_filename()+".png","bla",RST::target("/tmp/"+cs.get_filename()+".pdf")+RST::scale("200"));
					rst.text(out.get_header());
					rst.save(false,true);
					command.html_browser("/tmp/"+cs.get_filename()+".html");
					//if( my::get_yn("move the plot in the correct folder ?") ){
					//std::cerr<<"should implement this move"<<std::endl;
					//}
				}
			}
		}
	}
}
