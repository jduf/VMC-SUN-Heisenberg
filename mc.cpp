/*!@file mc.cpp */

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int i(0);

	Linux command;
	IOFiles* iof(NULL);
	System* sys;
	CreateSystem* cs(NULL);
	std::string fname;
	if(P.find("sim",i,false)){
		iof = new IOFiles(P.get<std::string>(i),false,false);
		Vector<double> tmp(*iof);
		sys = new System (*iof);
		sys->create_cluster(true);
		cs  = new CreateSystem(sys);
		cs->init(&tmp,NULL);
		fname = cs->get_filename();
	} else {
		sys = new System(P);
		cs  = new CreateSystem(sys);
		cs->init(NULL,&P);
		if(cs->get_status()==2 && P.find("obs",i,false)){
			Vector<unsigned int> obs(P.get_type(i)?P.get<std::vector<unsigned int> >(i):Vector<unsigned int>(1,P.get<unsigned int>(i)));
			for(unsigned int j(0);j<obs.size();j++){ cs->create_obs(obs(j)); }
		}
	}

	if(cs->get_status()==2){
		std::cout<<cs->get_system_info()<<std::endl;
		for(unsigned int j(0);j<cs->get_obs().size();j++){
			std::cout<<(j?"          ":"Compute : ")<<cs->get_obs()[j].get_name()<<std::endl;
		}
		std::cout<<std::endl;
		cs->create(false);
		if(cs->get_status()==1){
			if(!P.find("norun",i,false)){
				double dEoE(P.find("dEoE",i,false)?P.get<double>(i):0.001);
				unsigned int tmax(P.find("tmax",i,false)?P.get<unsigned int>(i):10);
				unsigned int nruns(P.find("nruns",i,false)?P.get<unsigned int>(i):omp_get_max_threads());
				unsigned int maxiter(P.find("maxiter",i,false)?P.get<unsigned int>(i):1);
				while(maxiter){
#pragma omp parallel for
					for(unsigned int j=0;j<nruns;j++){
#pragma omp critical
						std::cout<<"core number "<<j<<" "<<omp_get_thread_num()<<std::endl;
						MCSystem* mcsys(NULL);
						if(cs->use_complex()){
							if(cs->is_bosonic()){ mcsys = new SystemBosonic<std::complex<double> >(*dynamic_cast<const Bosonic<std::complex<double> >*>(cs->get_GenericSystem())); }
							else                { mcsys = new SystemFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs->get_GenericSystem())); }
						} else {
							if(cs->is_bosonic()){ mcsys = new SystemBosonic<double>(*dynamic_cast<const Bosonic<double>*>(cs->get_GenericSystem())); }
							else                { mcsys = new SystemFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs->get_GenericSystem())); }
						}
						if(!mcsys){ std::cerr<<__PRETTY_FUNCTION__<<" MCSystem was not constructed"<<std::endl; }
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

					maxiter--;
					if(std::abs(cs->get_obs()[0][0].get_dx()/cs->get_obs()[0][0].get_x())<dEoE){ maxiter=0; }
					else { std::cerr<<__PRETTY_FUNCTION__<<" : should rerun, bad precision : "<<cs->get_obs()[0][0].get_dx()/cs->get_obs()[0][0].get_x()<<" > dEoE = "<<dEoE<<std::endl; }
				}
				fname = Time().date("-") + "-" + cs->get_filename();
				command.mkpath(cs->get_path().c_str());
				if(iof){ delete iof; }
				iof = new IOFiles(cs->get_path() + fname +".jdbin",true,false);
				cs->save(*iof);
			} else {/* cs->check();*/ }
			cs->print(1);
		} else { std::cerr<<__PRETTY_FUNCTION__<<" : CreateSystem::create(&p,NULL) failed "<<std::endl; }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : CreateSystem::init(&p,NULL) failed (status="<<cs->get_status()<<")"<<std::endl; }
	if(cs->get_status()!=1){
		if(cs){ delete cs; cs = NULL; }
		std::cout<<RST::dash_line_<<std::endl;
	}

	if(P.find("d",i,false) && cs && iof){
		RSTFile rst("/tmp/",fname);
		IOSystem ios(fname,"","","","","/tmp/",&rst);
		cs->set_IOSystem(&ios);
		cs->display_results();

		rst.text(iof->get_header());
		rst.save(false,true);
		command(Linux::html_browser("/tmp/"+fname+".html"),true);
	}

	if(cs) { delete cs; }
	if(sys){ delete sys; }
	if(iof){ delete iof; }
}
