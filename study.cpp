/*! @file study.cpp */

#include "AnalyseMagnetization.hpp"
#include "AnalyseEnergy.hpp"

#include "AnalyseChain.hpp"
#include "AnalyseHoneycomb.hpp"

int main(int argc, char* argv[]){
	Linux command;
	command("[ -d sim ] && [ -d src ]");
	if(!command.status()){
		Parseur P(argc,argv);
		unsigned int i(0);
		std::string path(P.search("path",i)?P.get<std::string>(i):"");

		unsigned int what(P.get<unsigned int>("what"));
		
		if(!P.status() && what<4){
			switch(what){
				case 0: { AnalyseEnergy ana(path); }break;
				case 1: { AnalyseChain ana(path); }break;
				case 2: { AnalyseHoneycomb ana(path); }break;
				case 3: { AnalyseMagnetization ana(path); }break;
			}
		} else {
			std::cerr<<"study : unknown study"<<std::endl;
			std::cerr<<"study : unknown option what, options are :"<<std::endl; 
			std::cerr<<"         - Energy :        0"<<std::endl;
			std::cerr<<"         - Chain :         1"<<std::endl;
			std::cerr<<"         - Honeycomb     : 2"<<std::endl;
			std::cerr<<"         - Magnetization : 3"<<std::endl;
		}
	} else {
		std::cerr<<"study : must be run in a directory that contains 'sim' and 'src'"<<std::endl; 
	}
}
