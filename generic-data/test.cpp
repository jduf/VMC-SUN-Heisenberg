#include "Data.hpp"
#include "Matrix.hpp"
#include "Parseur.hpp"

int main(int argc,char* argv[]){
	Parseur P(argc,argv);
	Read r(P.get<std::string>("sim"));
	if(!P.status()){
		Input i;
		i.read<std::string>("wf",r);
		i.read<unsigned int>("N",r);
		i.read<unsigned int>("m",r);
		i.read<Matrix<unsigned int> >("sts",r);
		unsigned int n(10);
		i.set("nthreads",n);

		std::cout<<"N "<<i.get<unsigned int>("N")<<std::endl;
		std::cout<<"wf "<<i.get<std::string>("wf")<<std::endl;
		std::cout<<"sts "<<i.get<Matrix<unsigned int> >("sts")<<std::endl;
		std::cout<<"nthreads "<<i.get<unsigned int>("nthreads")<<std::endl;
		std::cout<<"rien "<<i.get<unsigned int>("fausse entrée")<<std::endl;
	}
}
