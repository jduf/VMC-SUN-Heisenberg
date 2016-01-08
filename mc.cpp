/*!  @file mc.cpp */

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"
#include <omp.h>

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int i(0);
	unsigned int tmax(P.get<unsigned int>("tmax"));
	unsigned int nruns(P.find("nruns",i,false)?P.get<unsigned int>(i):omp_get_max_threads());

	System* sys;
	CreateSystem* cs;
	if(P.find("sim",i,false)){
		IOFiles read(P.get<std::string>(i),false);
		Vector<double> tmp(read);
		sys = new System(read);
		cs  = new CreateSystem(sys);
		cs->init(&tmp,NULL);
	} else {
		if(!P.find("M",i,false)){
			std::vector<unsigned int> M(P.get<unsigned int>("N"),P.get<unsigned int>("n")*P.get<unsigned int>("m")/P.get<unsigned int>("N"));
			P.set("M",M);
		}
		sys = new System(P);
		cs  = new CreateSystem(sys);
		cs->init(NULL,&P);
		cs->set_obs(P.find("nobs",i,false)?P.get<int>(i):-1);
	}

	if(!P.locked()){
		if(cs->get_status()==2){
			cs->create(true);
			if(cs->get_status()==1){
#pragma omp parallel for
				for(unsigned int j=0;j<nruns;j++){
					MCSystem* mcsys(NULL);
					if(cs->use_complex()){
						if(cs->is_bosonic()){ mcsys = new SystemBosonic<std::complex<double> >(*dynamic_cast<const Bosonic<std::complex<double> >*>(cs->get_GenericSystem())); }
						else                { mcsys = new SystemFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs->get_GenericSystem())); }
					} else {
						if(cs->is_bosonic()){ mcsys = new SystemBosonic<double>(*dynamic_cast<const Bosonic<double>*>(cs->get_GenericSystem())); }
						else                { mcsys = new SystemFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs->get_GenericSystem())); }
					}
					if(!mcsys){ std::cout<<__PRETTY_FUNCTION__<<" MCSystem was not constructed"<<std::endl; }
					else {

						MonteCarlo sim(mcsys,tmax);
						sim.thermalize(1e7);
						sim.run(1e8);

#pragma omp critical(System__merge)
						{ sys->merge(mcsys); }

						delete mcsys;
					}
				}

				sys->complete_analysis(1e-5);
				std::cout<<sys->get_energy()<<std::endl;

				Linux command;
				command.mkpath(cs->get_path().c_str());
				IOFiles out(cs->get_path() + cs->get_filename()+".jdbin",true);
				cs->save_param(out);
				cs->save_input(out);
				sys->save_output(out);

				if(P.find("d",i,false)){
					CreateSystem tmp(sys);
					tmp.init(NULL,&P);
					tmp.set_obs(P.find("nobs",i,false)?P.get<int>(i):-1);
					RSTFile rst("/tmp/",tmp.get_filename());
					IOSystem ios(tmp.get_filename(),"","","","","/tmp/",&rst);
					tmp.set_IOSystem(&ios);

					tmp.display_results();

					rst.text(out.get_header());
					rst.save(false,true);
					command(Linux::html_browser("/tmp/"+cs->get_filename()+".html"),true);
					//if( my::get_yn("move the plot in the correct folder ?") ){
					//std::cerr<<"should implement this move"<<std::endl;
					//}
				}
			} else { std::cout<<__PRETTY_FUNCTION__<<" : CreateSystem::create(&p,NULL) failed "<<std::endl; }
		} else { std::cout<<__PRETTY_FUNCTION__<<" : CreateSystem::init(&p,NULL) failed "<<std::endl; }
	} else { std::cout<<__PRETTY_FUNCTION__<<" : Parseur locked"<<std::endl; }

	delete sys;
	delete cs;
}
