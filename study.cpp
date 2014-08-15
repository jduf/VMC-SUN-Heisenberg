/*! @file study.cpp */

#include "AnalyseMagnetization.hpp"
#include "AnalyseChain.hpp"
#include "AnalyseEnergy.hpp"

int main(int argc, char* argv[]){
	Linux command;
	command("[ -d sim ] && [ -d src ]");
	if(!command.status()){
		Parseur P(argc,argv);
		unsigned int i(0);
		std::string path("");
		Analyse* analyse(NULL);

		unsigned int what(P.get<unsigned int>("what"));
		if(P.search("path",i)){path = P.get<std::string>(i);}
		if(!P.status() && what<4){
			switch(what){
				case 0: { analyse = new AnalyseEnergy("sim/"); }break;
				case 1: { analyse = new AnalyseChain("sim/"); }break;
				case 2: { analyse = new AnalyseMagnetization("sim/"); }break;
			}
			analyse->go(path);
			delete analyse;
			analyse = NULL;
		} else {
			std::cerr<<"study : unknown study"<<std::endl;
			std::cerr<<"study : unknown option what, options are :"<<std::endl; 
			std::cerr<<"         - Energy :        0"<<std::endl;
			std::cerr<<"         - Chain :         1"<<std::endl;
			std::cerr<<"         - Magnetization : 2"<<std::endl;
		}
	} else {
		std::cerr<<"study : must be run in a directory that contains 'sim' and 'src'"<<std::endl; 
	}
}
