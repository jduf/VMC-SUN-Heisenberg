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
	Parseur P(argc,argv);
	std::string sim(P.check_get("sim",std::string("sim/")));
	unsigned int max_level(P.check_get("ml",10));
	/*!bash_file
	  = 0 : execute the linux commands on the fly
	  = 1 : save the linux commands in a bash file
	  = 2 : save the linux commands in a bash file and execute them when
	  Analyse is destructed*/
	unsigned int bash_file(P.check_get("bash_file",0));
	unsigned int display(P.check_get("display",0));
	switch(P.check_get("what",666)){
		case 0: { AnalyseEnergy ana(sim,max_level,bash_file); }break;
		case 1: { AnalyseMin ana(sim,max_level,bash_file); }break;
		case 2: { AnalyseExtract ana(sim,max_level,bash_file,display); }break;
		case 3: { AnalyseSystematic ana(sim,max_level,bash_file); }break;
		case 4: { AnalyseChain ana(sim,max_level,bash_file); }break;
		case 5: { AnalyseHoneycomb ana(sim,max_level,bash_file); }break;
		case 6: { AnalyseLadder ana(sim,max_level,bash_file); }break;
		case 7: { AnalyseSystematicLadder ana(sim,max_level,bash_file,display); }break;
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
}
