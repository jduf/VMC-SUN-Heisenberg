/*!@file check.cpp */

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"
#include "MCSim.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::cout<<"############# Init System #################"<<std::endl;
	unsigned int i(0);
	unsigned int tmax(P.find("tmax",i,false)?P.get<unsigned int>(i):10);
	System s(P);
	std::cout<<"############# Init CreateSystem ###########"<<std::endl;
	CreateSystem cs(&s);
	switch(P.find("what",i,true)?P.get<unsigned int>(i):666){
		case 0:/*call CreateSystem::init*/
			{
				std::cout<<"############# Init GenericSystem ##########"<<std::endl;
				cs.init(NULL,&P);
			} break;
		case 1:/*call CreateSystem::init and check*/
			{
				std::cout<<"############# Init GenericSystem ##########"<<std::endl;
				cs.init(NULL,&P);
				std::cout<<"############# Check #######################"<<std::endl;
				if(cs.get_status()==2){ cs.check(); }
			} break;
		case 2:/*call CreateSystem::init create and check*/
			{
				std::cout<<"############# Init GenericSystem ##########"<<std::endl;
				cs.init(NULL,&P);
				if(cs.get_status()==2){
					std::cout<<"############# Create GenericSystem ########"<<std::endl;
					cs.create();
					std::cout<<"############# Check #######################"<<std::endl;
					cs.check();
				}
			} break;
		case 3:/*call CreateSystem::init create, MonteCarlo::run*/
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
					MonteCarlo sim(S,tmax);
					sim.thermalize(1e6);
					std::cout<<"############# Run Monte Carlo #############"<<std::endl;
					sim.run();
					S->complete_analysis(1e-5);
					std::cout<<S->get_energy()<<std::endl;
					delete S;
				}
			} break;
		case 4:/*call CreateSystem::(init,create), MCSystem, MonteCarlo::run and save*/
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
					MonteCarlo sim(S,tmax);
					sim.thermalize(1e6);
					std::cout<<"############# Run Monte Carlo #############"<<std::endl;
					sim.run();
					S->complete_analysis(1e-5);
					std::cout<<S->get_energy()<<std::endl;
					std::cout<<"############# Save MCSystem ###############"<<std::endl;
					IOFiles out("check.jdbin",true);
					S->write(out);
					delete S;
				}
			} break;
		case 5:/*call load MCSystem, MonteCarlo::run and save*/
			{
				MCSystem* S(NULL);
				std::cout<<"############# Load MCSystem ###############"<<std::endl;
				IOFiles in("check.jdbin",false);
				if(cs.use_complex()){
					S = new SystemFermionic<std::complex<double> >(in);
				} else {
					S = new SystemFermionic<double>(in);
				}
				std::cout<<"############# Init Monte Carlo ############"<<std::endl;
				MonteCarlo sim(S,tmax);
				sim.thermalize(10);
				std::cout<<"############# Run Monte Carlo #############"<<std::endl;
				sim.run();
				S->complete_analysis(1e-5);
				std::cout<<S->get_energy()<<std::endl;
				std::cout<<"############# Save MCSystem ###############"<<std::endl;
				IOFiles out("check-1.jdbin",true);
				S->write(out);
				delete S;
			} break;
		case 6:/*use MCSim to run two sim then, merge them*/
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
				mcsim2->copy_S(mcsim);
				mcsim2->set_obs(cs.get_obs()[0]);
				mcsim2->run(1e6,4);
				mcsim2->complete_analysis(1e-5);

				std::cout<<mcsim->get_energy()<<std::endl;
				std::cout<<mcsim2->get_energy()<<std::endl;
				mcsim->merge(mcsim2);
				mcsim->complete_analysis(1e-5);
				std::cout<<mcsim->get_energy()<<std::endl;
			} break;
		case 7:/*check symmetries*/
			{
				Vector<double> t_ref(P.get<std::vector<double> >("t"));
				Vector<double> mu_ref(P.get<std::vector<double> >("mu"));
				Vector<double> param_ref(t_ref.size()+mu_ref.size());
				for(unsigned int i(0);i<t_ref.size();i++){ param_ref(i) = t_ref(i); }
				for(unsigned int i(0);i<mu_ref.size();i++){ param_ref(t_ref.size()+i) = mu_ref(i); }
				cs.init(NULL,&P);
				std::vector<Matrix<int> > all_sym;
				cs.get_wf_symmetries(all_sym);
				for(unsigned int j(0);j<all_sym.size();j++){
					Vector<double> param(param_ref);
					Matrix<int> sym(all_sym[j]);
					for(unsigned int i(0);i<sym.row();i++){
						if(sym(i,1)<0){ param(sym(i,0)) = sym(i,2)*66; }
						else { param(sym(i,0)) = sym(i,2)*param(sym(i,1)); }
					}
					for(unsigned int i(0);i<param.size();i++){
						if( my::are_equal(std::abs(param(i)),66) ){
							std::cout<<(param(i)>=0?" A ":"-A ");
						} else {
							std::cout<<(param(i)>=0?" ":"")<<param(i)<<" ";
						}
					}
					std::cout<<std::endl;
				}
			} break;
		default:
			{
				std::cerr<<__PRETTY_FUNCTION__<<" : unknown option 'what', options are :"<<std::endl;
				std::cerr<<"    - init                        : 0"<<std::endl;
				std::cerr<<"    - init + check                : 1"<<std::endl;
				std::cerr<<"    - init + create + check       : 2"<<std::endl;
				std::cerr<<"    - init + create + run         : 3"<<std::endl;
				std::cerr<<"    - init + create + run + write : 4"<<std::endl;
				std::cerr<<"    - load + run + rewrite        : 5"<<std::endl;
				std::cerr<<"    - MCSim                       : 6"<<std::endl;
				std::cerr<<"    - load + check_symmetries     : 7"<<std::endl;
			}
	}
}
