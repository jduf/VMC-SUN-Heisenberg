/*!  @file mc.cpp */

#include "Read.hpp"
#include "Parseur.hpp"
#include "Gnuplot.hpp"
#include "Directory.hpp"
#include "Linux.hpp"

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::string directory_name(P.get<std::string>("0"));
	Directory d;
	Linux command;
	if(directory_name == "."){ directory_name = command.pwd(); }
	else { directory_name = command.pwd()+directory_name; }
	d.search_file_ext(".jdbin",directory_name,false);
	Write runs("runs.dat");
	Write mean("mean.dat");
	for(unsigned int i(0); i<d.size();i++){
		Read r(d[i]);
		unsigned int code;
		r>>code;
		switch(code){
			case 0:{
					   unsigned int nruns;
					   unsigned int N;
					   unsigned int m;
					   unsigned int n;
					   int bc;
					   double param;
					   double E;
					   double DeltaE;
					   unsigned int Nsteps;
					   unsigned int status;
					   Vector<double> corr;

					   r>>nruns>>N>>m>>n>>bc>>param;
					   for(unsigned int i(0);i<nruns;i++){
						   r>>E>>DeltaE>>Nsteps>>status>>corr;
						   runs<<param<<" "<<E<<" "<<DeltaE<<Write::endl;
					   }
				   }break;
			case 1:{
					   unsigned int nruns;
					   unsigned int N;
					   unsigned int m;
					   unsigned int n;
					   int bc;
					   double param;
					   double E;
					   double DeltaE;

					   r>>nruns>>N>>m>>n>>bc>>param;
					   do{
						   r>>param>>E>>DeltaE;
						   mean<<param<<" "<<E<<" "<<DeltaE<<Write::endl;
					   } while ( !r.eof() );
				   }break;
			default:{std::cerr<<"The program that created the file is unknown"<<std::endl;}
		}

		Gnuplot gp("plot","plot");
		gp.add_plot_param("'runs.dat' u 1:2:3 w errorbars, ");
		gp.add_plot_param("'trash.dat' u 1:2, ");
		gp.add_plot_param("'mean.dat' u 1:2:3 w errorbars ");
	}
	std::cerr<<"modifer Read->eof et peut être tout ce qui va avec"<<std::endl;
}
