/*!  @file study.cpp */

#include "AnalyseMagnetization.hpp"
#include "AnalyseParameter.hpp"

int main(int argc, char* argv[]){
	Linux command;
	command("[ -d sim ] && [ -d src ]");
	if(!command.status()){
		Parseur P(argc,argv);
		unsigned int what(P.get<unsigned int>("what"));
		unsigned int i(0);
		switch(what){
			case 1:
				{
					AnalyseParameter analyse("sim/"); 
					if(!P.search("dir",i)){analyse.go("");}
					else {analyse.go(P.get<std::string>(i));}
				}break;
			case 2:
				{
					AnalyseMagnetization analyse("sim/"); 
					if(!P.search("dir",i)){analyse.go("");}
					else {analyse.go(P.get<std::string>(i));}
				}break;
			default:
				{
					std::cerr<<"study : unknown study"<<std::endl;
				}break;
		}
	} else {
		std::cerr<<"study : must be run in a directory that contains 'sim' and 'src'"<<std::endl;
	}
}
