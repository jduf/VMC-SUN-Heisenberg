/*!@file mcbi.cpp */

#include "BiSystem.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int i(0);
	unsigned int tmax(P.get<unsigned int>("tmax"));
	unsigned int nruns(P.find("nruns",i,false)?P.get<unsigned int>(i):omp_get_max_threads());

	if(P.find("load",i,false)){
		IOFiles r(P.get<std::string>(i),false);
		BiSystem bs(r);

		bs.compute_dE();
		std::cout<<bs.get_E()<<std::endl;
		std::cout<<bs.get_dE()<<std::endl;
	} else {
		if(!P.find("M",i,false)){
			std::vector<unsigned int> M(P.get<unsigned int>("N"),P.get<unsigned int>("n")*P.get<unsigned int>("m")/P.get<unsigned int>("N"));
			P.set("M",M);
		}

		BiSystem bs(P);
		IOFiles r(P.get<std::string>("params"),false);
		Matrix<double> MP(P.get<unsigned int>("nwfs"),4);
		if(!P.locked()){
			r>>MP;
			Vector<double> param(MP.col());
			for(unsigned int i(0);i<MP.row();i++){
				for(unsigned int j(0);j<MP.col();j++){ param(j) = MP(i,j); }
				bs.add_new_param(param);
			}
			bs.run(nruns,tmax);
			bs.compute_E();
			std::cout<<bs.get_E()<<std::endl;
			std::cout<<bs.get_H()<<std::endl;
			std::cout<<bs.get_O()<<std::endl;
			bs.save();
		} else { std::cout<<__PRETTY_FUNCTION__<<" : Parseur locked"<<std::endl; }
	}
}
