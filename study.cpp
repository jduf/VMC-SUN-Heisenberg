/*!  @file study.cpp */

#include "AnalyseMagnetization.hpp"
#include "AnalyseParameter.hpp"

int main(int argc, char* argv[]){
	//if(argc==1){ AnalyseMagnetization a(""); }
	//else { AnalyseMagnetization a(argv[1]); }
	
	Parseur P(argc,argv);
	unsigned int what(P.get<unsigned int>("what"));
	unsigned int i(0);
	switch(what){
		case 1:
			{
				AnalyseParameter analyse; 
				if(P.search("dir",i)){analyse.go("");}
				else {analyse.go(P.get<std::string>(i));}
			}
		case 2:
			{
				AnalyseMagnetization analyse; 
				if(P.search("dir",i)){analyse.go("");}
				else {analyse.go(P.get<std::string>(i));}
			}
	}

}
