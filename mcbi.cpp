/*!@file mcbi.cpp */

#include "BiSystem.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int i(0);

	if(P.find("load",i,false)){
		IOFiles r(P.get<std::string>(i),false);
		BiSystem bs(r);

		//bs.compute_dE();
		//std::cout<<bs.get_E()<<std::endl;
		//std::cout<<bs.get_dE()<<std::endl;
		//std::cout<<std::endl;
		//std::cout<<bs.get_H()<<std::endl;
		//std::cout<<bs.get_dH()<<std::endl;
		//std::cout<<std::endl;
		//std::cout<<bs.get_O()<<std::endl;
		//std::cout<<bs.get_dO()<<std::endl;
		//bs.run(nruns,tmax);
		//bs.compute_E();
		//bs.save();
		bs.study();
	} else {
		unsigned int tmax(P.get<unsigned int>("tmax"));
		unsigned int nruns(P.find("nruns",i,false)?P.get<unsigned int>(i):omp_get_max_threads());

		BiSystem bs(P);
		if(!P.locked()){
			if(P.find("chain",i,false)){
				double t(P.get<double>("t"));
				for(unsigned int s(0);s<4;s++){
					Vector<double> param(4,1);
					param(s) = t;
					bs.add_new_param(param);
				}
			}
			if(P.find("ladder",i,false)){
				Vector<double> param(20,0);
				double t(P.get<double>("t"));
				for(unsigned int s(0);s<4;s++){
					for(unsigned int i(0);i<12;i++){
						switch(i%3){
							case 0: { param(i) = 1.00; } break;
							case 1: { param(i) = 0.15; } break;
							case 2: { param(i) =-1.00; } break;
						}
					}
					param(s*3)  = t;
					param(s*3+2)=-t;
					bs.add_new_param(param);
				}
			}
			if(P.find("honeycomb",i,false)){
				double t(P.get<double>("td"));
				Vector<double> param(2,0);
				param(0) = t;
				for(unsigned int fc(0);fc<3;fc++){
					param(1) = fc;
					bs.add_new_param(param);
				}
			}
			if(P.find("boundary",i,false)){
				//Vector<double> param(8,0);
				//param(0) = 1;
				//param(1) = 1;
				//param(2) = 1;
				//param(3) = 0.99999999;
				//bs.add_new_param(param);
				//param(3) =-0.99999999;
				//bs.add_new_param(param);
				Vector<double> param(20,0);
				param(0) = 1;
				param(1) = 0.15;
				param(2) =-1;
				param(3) = 1;
				param(4) = 0.15;
				param(5) =-1;
				param(6) = 1;
				param(7) = 0.15;
				param(8) =-1;
				param(9) = 1;
				param(10)= 0.15;
				param(11)=-1;
				param(12)= 0.9999999;
				bs.add_new_param(param);
				param(12)=-0.9999999;
				bs.add_new_param(param);
			}

			if(P.find("square",i,false)){
				for(unsigned int s(0);s<5;s++){
					Vector<double> mu(5,0);
					mu(s) = P.get<double>("mu");
					bs.add_new_param(mu);
				}
			}
			bs.run(nruns,tmax);
			bs.compute_E();
			std::cout<<bs.get_E()<<std::endl;
			std::cout<<std::endl;
			std::cout<<bs.get_H()<<std::endl;
			std::cout<<std::endl;
			std::cout<<bs.get_O()<<std::endl;
			bs.save();
		} else { std::cout<<__PRETTY_FUNCTION__<<" : Parseur locked"<<std::endl; }
	}
}
