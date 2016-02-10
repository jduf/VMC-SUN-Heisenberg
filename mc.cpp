/*!@file mc.cpp */

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"

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
		sys->print(1);
		cs  = new CreateSystem(sys);
		cs->init(&tmp,NULL);
	} else {
		sys = new System(P);
		cs  = new CreateSystem(sys);
		cs->init(NULL,&P);
		cs->set_obs(P.find("nobs",i,false)?P.get<int>(i):-1);
	}

	if(!P.locked()){
		std::cout<<cs->get_system_info()<<std::endl;
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
						sim.thermalize(1e6);
						sim.run();

#pragma omp critical(System__merge)
						{ cs->merge(mcsys); }

						delete mcsys;
					}
				}
				cs->complete_analysis(1e-5);
				cs->print(2);

				Linux command;
				command.mkpath(cs->get_path().c_str());
				std::string fname(cs->get_filename() + "-" + Time().date("-") );
				IOFiles out(cs->get_path() + fname +".jdbin",true);
				cs->save(out);

				if(P.find("d",i,false)){
					RSTFile rst("/tmp/",fname);
					IOSystem ios(fname,"","","","","/tmp/",&rst);
					cs->set_IOSystem(&ios);
					cs->display_results();

					rst.text(out.get_header());
					rst.save(false,true);
					command(Linux::html_browser("/tmp/"+fname+".html"),true);
					//if( my::get_yn("move the plot in the correct folder ?") ){
					//std::cerr<<"should implement this move"<<std::endl;
					//}
				}
			} else { std::cout<<__PRETTY_FUNCTION__<<" : CreateSystem::create(&p,NULL) failed "<<std::endl; }
		} else { std::cout<<__PRETTY_FUNCTION__<<" : CreateSystem::init(&p,NULL) failed "<<std::endl; }
		if(cs->get_status()>1){ std::cout<<std::string(35,'-')<<std::endl; }
	} else { std::cout<<__PRETTY_FUNCTION__<<" : Parseur locked"<<std::endl; }

	delete cs;
	delete sys;
}
