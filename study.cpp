/*! @file study.cpp */

#include "AnalyseEnergy.hpp"
#include "AnalyseMin.hpp"
#include "AnalyseExtract.hpp"
#include "AnalyseSystematic.hpp"
#include "AnalyseSystematicCoupling.hpp"
#include "AnalyseChain.hpp"
#include "AnalyseHoneycomb.hpp"
#include "AnalyseLadder.hpp"
#include "AnalyseLadderPSO.hpp"

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
		case 0: { AnalyseEnergy(sim,max_level,bash_file,P.get<unsigned int>("ref")); }break;
		case 1: { AnalyseMin(sim,max_level,bash_file,P.get<unsigned int>("ref")); }break;
		case 2: { AnalyseExtract(sim,max_level,bash_file,display,P.get<unsigned int>("ref")); }break;
		case 3: { AnalyseSystematic(sim,max_level,bash_file,P.get<unsigned int>("ref")); }break;
		case 4: { AnalyseSystematicCoupling(sim,max_level,bash_file,display,P.get<unsigned int>("ref")); }break;
		case 5: { AnalyseChain(sim,max_level,bash_file); }break;
		case 6: { AnalyseHoneycomb(sim,max_level,bash_file); }break;
		case 7: { AnalyseLadder(sim,max_level,bash_file); }break;
		case 8: { AnalyseLadderPSO(sim,max_level,bash_file,display); }break;
		default:{
					std::cerr<<__PRETTY_FUNCTION__<<" : unknown option 'what', options are :"<<std::endl;
					std::cerr<<"  - Energy            : 0"<<std::endl;
					std::cerr<<"  - Min               : 1"<<std::endl;
					std::cerr<<"  - Extract           : 2"<<std::endl;
					std::cerr<<"  - Systematic        : 3"<<std::endl;
					std::cerr<<"  - SystematicCoupling: 4"<<std::endl;
					std::cerr<<"  - Chain             : 5"<<std::endl;
					std::cerr<<"  - Honeycomb         : 6"<<std::endl;
					std::cerr<<"  - Ladder            : 7"<<std::endl;
					std::cerr<<"  - LadderPSO         : 9"<<std::endl;
				}
	}
}
