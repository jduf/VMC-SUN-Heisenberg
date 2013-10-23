#include "Parseur.hpp"
#include "Read.hpp"
#include "PSTricks.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::string filename(P.get<std::string>("0"));
	if(!P.status()){
		Read file(filename);
		std::string s("");
		file>>s;

		std::size_t a(0);
		std::size_t b(0);

		std::string begin("\\begin{pspicture}");
		std::string end("\\end{pspicture}");

		while(a != std::string::npos){
			a=s.find(begin,b);
			if(a != std::string::npos){
				b=s.find(end,a)+end.size();

				std::string nameline(s.substr(a,s.find('\n',a)-a));
				std::string outname(nameline.substr(nameline.find('%')+1));

				std::cout<<outname<<std::endl;

				if(outname!=nameline){
					PSTricks ps(outname);
					ps.add(s.substr(a,b-a));
				};
			}
		}
	} else {
		std::cerr<<"pstricks2pdf : You need to give a .tex file"<<std::endl;
	}
}
