/*!@file check.cpp */

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"
#include "MCSim.hpp"

#include "AnalyseExtract.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int what(P.check_get("what",0));
	if(what<6){
		std::cout<<"############# Init System #################"<<std::endl;
		System s(P);
		if(s.get_status() == 4){
			std::cout<<"############# Init CreateSystem ###########"<<std::endl;
			CreateSystem cs(&s);
			switch(what){
				case 0:/*call CreateSystem::init and check*/
					{
						std::cout<<"############# Init GenericSystem ##########"<<std::endl;
						cs.init(NULL,&P);
						std::cout<<"############# Check #######################"<<std::endl;
						if(cs.get_status()==2){
							std::cout<<cs.get_system_info()<<std::endl;
							cs.check();
						}
					}break;
				case 1:/*call CreateSystem::init and create*/
					{
						std::cout<<"############# Init GenericSystem ##########"<<std::endl;
						cs.init(NULL,&P);
						if(cs.get_status()==2){
							std::cout<<"############# Create GenericSystem ########"<<std::endl;
							cs.create();
							if(cs.get_status()!=1){ s.print(true); }
						}
					}break;
				case 2:/*call CreateSystem::init create, MonteCarlo::run*/
					{
						std::cout<<"############# Init GenericSystem ##########"<<std::endl;
						cs.init(NULL,&P);
						cs.create_obs(0);
						if(cs.get_status()==2){
							std::cout<<"############# Create GenericSystem ########"<<std::endl;
							cs.create();
							std::cout<<"############# Create MCSystem #############"<<std::endl;
							MCSystem* S(NULL);
							if(cs.use_complex()){
								S = new SystemFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs.get_GenericSystem()));
							} else {
								S = new SystemFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs.get_GenericSystem()));
							}
							std::cout<<"############# Init Monte Carlo ############"<<std::endl;
							MonteCarlo sim(S,P.check_get("tmax",10));
							sim.thermalize(1e6);
							std::cout<<"############# Run Monte Carlo #############"<<std::endl;
							sim.run();
							S->complete_analysis(1e-5);
							std::cout<<S->get_energy()<<std::endl;
							delete S;
						}
					}break;
				case 3:/*call CreateSystem::(init,create), MCSystem, MonteCarlo::run and save*/
					{
						std::cout<<"############# Init GenericSystem ##########"<<std::endl;
						cs.init(NULL,&P);
						cs.create_obs(0);
						if(cs.get_status()==2){
							std::cout<<"############# Create GenericSystem ########"<<std::endl;
							cs.create();
							std::cout<<"############# Create MCSystem #############"<<std::endl;
							MCSystem* S(NULL);
							if(cs.use_complex()){
								S = new SystemFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs.get_GenericSystem()));
							} else {
								S = new SystemFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs.get_GenericSystem()));
							}
							S->set_obs(cs.get_obs()[0]);
							std::cout<<"############# Init Monte Carlo ############"<<std::endl;
							MonteCarlo sim(S,P.check_get("tmax",10));
							sim.thermalize(1e6);
							std::cout<<"############# Run Monte Carlo #############"<<std::endl;
							sim.run();
							S->complete_analysis(1e-5);
							std::cout<<S->get_energy()<<std::endl;
							std::cout<<"############# Save MCSystem ###############"<<std::endl;
							IOFiles out("check.jdbin",true,false);
							S->write(out);
							delete S;
						}
					}break;
				case 4:/*call load MCSystem, MonteCarlo::run and save*/
					{
						MCSystem* S(NULL);
						std::cout<<"############# Load MCSystem ###############"<<std::endl;
						IOFiles in("check.jdbin",false,false);
						if(cs.use_complex()){
							S = new SystemFermionic<std::complex<double> >(in);
						} else {
							S = new SystemFermionic<double>(in);
						}
						std::cout<<"############# Init Monte Carlo ############"<<std::endl;
						MonteCarlo sim(S,P.check_get("tmax",10));
						sim.thermalize(10);
						std::cout<<"############# Run Monte Carlo #############"<<std::endl;
						sim.run();
						S->complete_analysis(1e-5);
						std::cout<<S->get_energy()<<std::endl;
						std::cout<<"############# Save MCSystem ###############"<<std::endl;
						IOFiles out("check-1.jdbin",true,false);
						S->write(out);
						delete S;
					}break;
				case 5:/*use MCSim to run two sim then, merge them*/
					{
						cs.init(NULL,&P);
						cs.create();
						cs.create_obs(0);

						std::shared_ptr<MCSim> mcsim(std::make_shared<MCSim>(P.get<std::vector<double> >("t")));
						System s(P);
						mcsim->create_S(&s);
						mcsim->set_obs(cs.get_obs()[0]);

						mcsim->run(1e6,2);
						mcsim->complete_analysis(1e-5);

						std::shared_ptr<MCSim> mcsim2(std::make_shared<MCSim>(P.get<std::vector<double> >("t")));
						mcsim2->copy_clear_S(mcsim);
						mcsim2->set_obs(cs.get_obs()[0]);
						mcsim2->run(1e6,4);
						mcsim2->complete_analysis(1e-5);

						std::cout<<mcsim->get_energy()<<std::endl;
						std::cout<<mcsim2->get_energy()<<std::endl;
						mcsim->merge(mcsim2);
						mcsim->complete_analysis(1e-5);
						std::cout<<mcsim->get_energy()<<std::endl;
					}break;
				default:
					{
						std::cerr<<__PRETTY_FUNCTION__<<" : unknown option 'what', options are :"<<std::endl;
						std::cerr<<"    - init + check                : 0"<<std::endl;
						std::cerr<<"    - init + create               : 1"<<std::endl;
						std::cerr<<"    - init + create + run         : 2"<<std::endl;
						std::cerr<<"    - init + create + run + write : 3"<<std::endl;
						std::cerr<<"    - load + run + rewrite        : 4"<<std::endl;
						std::cerr<<"    - MCSim                       : 5"<<std::endl;
					}
			}
		}
	} else {/*extract a full minimisation simulation */

		//IOFiles read(P.get<std::string>("sim"),false,false);

		//VMCExtract min(read,1e3,1e4);
		//List<MCSim> kept_samples;
		//List<MCSim>::Node* target(min.analyse("./","test",kept_samples));
		//Interpolation<Vector<double> > inter(12);
		//if(!kept_samples.size()){ std::cerr<<__PRETTY_FUNCTION__<<" : need to have at least one sample in 'kept_samples_'"<<std::endl; }
		//else {
			//kept_samples.set_target();
			//unsigned int l(target->get()->get_param().size());
			//while(kept_samples.target_next()){
				////if((target->get()->get_param()-kept_samples.get().get_param()).norm()/l<0.1){
					////std::cout<<kept_samples.get().get_param()<<" "<<kept_samples.get().get_energy()<<std::endl;
				////}
			//}
		//}
	}
}
