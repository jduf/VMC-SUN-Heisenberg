#include "Linux.hpp" 
#include "Time.hpp" 
#include "Parseur.hpp" 
#include "Directory.hpp" 
#include "RSTFile.hpp" 

#include <omp.h>

int main(int argc,char* argv[]){
	unsigned int min_args(2);
	Parseur P(argc,argv);
	std::string ext(P.get<std::string>("e"));

	unsigned int max_proc(P.get<unsigned int>("t"));

	if(!P.status()){
		Directory D;
		Linux L;
		std::string save_in("sim-info/");
		L("mkdir " + save_in);

		D.search_file_ext(ext,L.pwd(),false,false);
		D.sort();

		std::string command("");
		for(unsigned int i(min_args); i<P.size();i++){
			command += P.get<std::string>(i) + " ";
		}

		Time t;
		std::string fname("");
		fname = save_in + my::tostring(t.year()) 
			+ "-" + my::tostring(t.month()) 
			+ "-" + my::tostring(t.day())
			+ "_" + my::tostring(t.hour()) + my::tostring(t.min());
#pragma omp parallel for num_threads(max_proc) 
		for(unsigned int i=0;i<D.size();i++) {
			std::cout<<command<<D.get_name(i)<<D.get_ext(i)<<std::endl;
			L(command + D.get_name(i) + D.get_ext(i) + "> " + D.get_name(i) + "_" + fname + "_" + my::tostring(i)+".log 2> " + D.get_name(i) + "_" + fname + "_" + my::tostring(i)+".err.log");
		}

		RSTFile rst("./","SIM");
		rst.title("bla",'=');

		L("firefox SIM.html");
	}
}

