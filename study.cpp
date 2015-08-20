/*! @file study.cpp */

#include "AnalyseMagnetization.hpp"
#include "AnalyseEnergy.hpp"

#include "AnalyseChain.hpp"
#include "AnalyseLadder.hpp"
#include "AnalyseHoneycomb.hpp"

int main(int argc, char* argv[]){
	Linux command;
	command("[ -d sim ] && [ -d src ]",false);
	if(!command.status()){
		Parseur P(argc,argv);
		unsigned int i(0);
		std::string path(P.find("path",i,false)?P.get<std::string>(i):"");
		unsigned int max_level(P.find("ml",i,false)?P.get<unsigned int>(i):10);
		switch(P.find("what",i,true)?P.get<unsigned int>(i):666){
			case 0: { AnalyseEnergy ana(path,max_level); }break;
			case 1: { AnalyseChain ana(path,max_level); }break;
			case 2: { AnalyseLadder ana(path,max_level); }break;
			case 3: { AnalyseHoneycomb ana(path,max_level); }break;
			case 4: { AnalyseMagnetization ana(path,max_level); }break;
			default:{
						std::cerr<<"study : unknown option 'what', options are :"<<std::endl;
						std::cerr<<"      - Energy        : 0"<<std::endl;
						std::cerr<<"      - Chain         : 1"<<std::endl;
						std::cerr<<"      - Ladder        : 2"<<std::endl;
						std::cerr<<"      - Honeycomb     : 3"<<std::endl;
						std::cerr<<"      - Magnetization : 4"<<std::endl;
					}
		}
	} else {
		std::cerr<<"study : must be run in a directory that contains 'sim' and 'src'"<<std::endl; 
	}
}
