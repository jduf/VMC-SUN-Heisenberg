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
	if(P.find("sim",i)){
		iof = new IOFiles(P.get<std::string>(i),false,false);
		Vector<double> tmp(*iof);
		sys = new System (*iof);
		sys->create_cluster(true);
		cs  = new CreateSystem(sys);
		cs->init(&tmp,NULL);
		fname = cs->get_filename();
	} else {
		sys = new System(P);
		if(sys->get_status() == 4){
			cs  = new CreateSystem(sys);
			cs->init(NULL,&P);
			if(cs->get_status()==2 && P.find("obs",i)){
				Vector<unsigned int> obs(P.get_type(i)?P.get<std::vector<unsigned int> >(i):Vector<unsigned int>(1,P.get<unsigned int>(i)));
				for(unsigned int j(0);j<obs.size();j++){ cs->create_obs(obs(j)); }
			}
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
			if(!P.find("norun")){
				double dEoE(P.check_get("dEoE",1e-5));
				unsigned int tmax(P.check_get("tmax",10));
				unsigned int nruns(P.check_get("nruns",omp_get_max_threads()));
				unsigned int maxiter(P.check_get("maxiter",1));
				while(maxiter){
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
					cs->complete_analysis();

					maxiter--;
					if(std::abs(cs->get_dEoE())<dEoE){ maxiter=0; }
					else { std::cerr<<__PRETTY_FUNCTION__<<" : should rerun, bad precision : "<<cs->get_dEoE()<<" > dEoE = "<<dEoE<<std::endl; }
				}
				fname = Time().date("-") + "-" + cs->get_filename();
				command.mkpath(cs->get_path().c_str());
				if(iof){ delete iof; }
				iof = new IOFiles(cs->get_path() + fname +".jdbin",true,false);
				cs->save(*iof);
			} else {/* cs->check();*/ }
			cs->print(true);
		} else { std::cerr<<__PRETTY_FUNCTION__<<" : CreateSystem::create(&p,NULL) failed "<<std::endl; }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : CreateSystem::init(&p,NULL) failed (status="<<cs->get_status()<<")"<<std::endl; }
	if(cs->get_status()!=1){
		if(cs){ delete cs; cs = NULL; }
		std::cout<<RST::dash_line_<<std::endl;
	}

	/*Display the results in a html browser*/
	if(P.find("d") && cs && iof){
		RSTFile rst("/tmp/",fname);
		IOSystem ios(fname,"","","","","/tmp/",&rst,false);
		cs->set_IOSystem(&ios);
		cs->display_results();

		rst.text(iof->get_header());
		rst.save(true,false,true);
		command(Linux::html_browser("/tmp/"+fname+".html"),true);
	}

	/*Print the results in a pdf and html files*/
	if(P.find("p") && cs && iof){
		RSTFile rst("./",fname+"-print");
		IOSystem ios(fname,"./","./","./","./","./",&rst,false);
		cs->set_IOSystem(&ios);
		cs->display_results();

		rst.text(iof->get_header());
		rst.save(true,true,true);

		command("rm " + fname + "*-print.tex ",true);
		command.mkdir(fname.c_str());
		command("mv " + fname + "* " + fname,true);
	}

	if(cs) { delete cs; }
	if(sys){ delete sys; }
	if(iof){ delete iof; }
}
