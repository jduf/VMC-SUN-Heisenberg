/*! @file study.cpp */

#include "AnalyseEnergy.hpp"

#include "AnalyseChain.hpp"
#include "AnalyseHoneycomb.hpp"
#include "AnalyseLadder.hpp"

#include "AnalyseMin.hpp"
#include "AnalyseExtract.hpp"
#include "AnalyseSystematic.hpp"
#include "AnalyseSystematicLadder.hpp"

int main(int argc, char* argv[]){
	Linux command;
	if(!command.status()){
		Parseur P(argc,argv);
		std::string path(P.check_get("path",std::string("")));
		std::string sim(P.check_get("sim",std::string("sim/")));
		unsigned int max_level(P.check_get("ml",10));
		unsigned int run_cmd(P.check_get("run",0));
		unsigned int display(P.check_get("display",0));
		switch(P.check_get("what",666)){
			case 0: { AnalyseEnergy ana(sim,path,max_level,run_cmd); }break;
			case 1: { AnalyseMin ana(sim,path,max_level,run_cmd); }break;
			case 2: { AnalyseExtract ana(sim,path,max_level,run_cmd,display); }break;
			case 3: { AnalyseSystematic ana(sim,path,max_level,run_cmd); }break;
			case 4: { AnalyseChain ana(sim,path,max_level,run_cmd); }break;
			case 5: { AnalyseHoneycomb ana(sim,path,max_level,run_cmd); }break;
			case 6: { AnalyseLadder ana(sim,path,max_level,run_cmd); }break;
			case 7: { AnalyseSystematicLadder ana(sim,path,max_level,run_cmd,display); }break;
			default:{
						std::cerr<<__PRETTY_FUNCTION__<<" : unknown option 'what', options are :"<<std::endl;
						std::cerr<<"   - Energy           : 0"<<std::endl;
						std::cerr<<"   - Min              : 1"<<std::endl;
						std::cerr<<"   - Extract          : 2"<<std::endl;
						std::cerr<<"   - Systematic       : 3"<<std::endl;
						std::cerr<<"   - Chain            : 4"<<std::endl;
						std::cerr<<"   - Honeycomb        : 5"<<std::endl;
						std::cerr<<"   - Ladder           : 6"<<std::endl;
						std::cerr<<"   - SystematicLadder : 7"<<std::endl;
					}
		}
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : must be run in a directory that contains 'sim' and 'src'"<<std::endl;
	}
}
