/*! @file study.cpp */

#include "AnalyseMagnetization.hpp"
#include "AnalyseEnergy.hpp"

#include "AnalyseChain.hpp"
#include "AnalyseHoneycomb.hpp"

#include "AnalyseMin.hpp"
#include "AnalyseExtract.hpp"

int main(int argc, char* argv[]){
	Linux command;
	if(!command.status()){
		Parseur P(argc,argv);
		unsigned int i(0);
		std::string path(P.find("path",i,false)?P.get<std::string>(i):"");
		std::string sim(P.find("sim",i,false)?P.get<std::string>(i):"sim/");
		unsigned int max_level(P.find("ml",i,false)?P.get<unsigned int>(i):10);
		unsigned int run_cmd(P.find("run",i,false)?P.get<unsigned int>(i):0);
		switch(P.find("what",i,true)?P.get<unsigned int>(i):666){
			case 0: { AnalyseEnergy ana(sim,path,max_level,run_cmd); }break;
			case 1: { AnalyseChain ana(sim,path,max_level,run_cmd); }break;
			case 2: { AnalyseHoneycomb ana(sim,path,max_level,run_cmd); }break;
			case 3: { AnalyseMagnetization ana(sim,path,max_level,run_cmd); }break;
			case 4: { AnalyseMin ana(sim,path,max_level,run_cmd); }break;
			case 5: { AnalyseExtract ana(sim,path,max_level,run_cmd,(P.find("dr",i,false)?P.get<bool>(i):false)); }break;
			default:{
						std::cerr<<__PRETTY_FUNCTION__<<" : unknown option 'what', options are :"<<std::endl;
						std::cerr<<"    - Energy        : 0"<<std::endl;
						std::cerr<<"    - Chain         : 1"<<std::endl;
						std::cerr<<"    - Honeycomb     : 2"<<std::endl;
						std::cerr<<"    - Magnetization : 3"<<std::endl;
						std::cerr<<"    - Min           : 4"<<std::endl;
					}
		}
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : must be run in a directory that contains 'sim' and 'src'"<<std::endl;
	}
}
