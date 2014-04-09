#include "Container.hpp"
#include "Matrix.hpp"
#include "Parseur.hpp"

int main(int argc,char* argv[]){
	Parseur P(argc,argv);
	FileParser file(P.get<std::string>("sim"));
	if(!P.status()){
		Container input;
		std::string wf;
		file.extract(wf);
		file.extract<unsigned int>("N",input);
		file.extract<unsigned int>("m",input);
		input.set("a",21.23);
		unsigned int n(10);
		input.set("nthreads",n);
		Matrix<double> mat(2,3,1.2);
		mat(1,2) = 5.2;
		input.set("mat",mat);

		std::cout<<"will output input"<<std::endl;
		//std::cout<<input<<std::endl;
		std::cout<<"N "<<input.get<unsigned int>("N")<<std::endl;
		std::cout<<"wf "<<wf<<std::endl;
		std::cout<<"nthreads "<<input.get<unsigned int>("nthreads")<<std::endl;
		//Container output(true);
		//output.set("N",input.get<unsigned int>("N"));
		//output.set("m",input.get<unsigned int>("m"));
		//for(unsigned int i(0);i<output.size();i++){
			//std::cout<<output.name(i)<<" ";
		//}
		//std::cout<<std::endl;
		//std::cout<<output<<std::endl;

		std::cout<<"copy input"<<std::endl;
		Container copyinput(input);
		std::cout<<"N "<<copyinput.get<unsigned int>("N")<<std::endl;
		std::cout<<"m "<<copyinput.get<unsigned int>("m")<<std::endl;
		std::cout<<"a "<<copyinput.get<double>("a")<<std::endl;
		std::cout<<"mat "<<std::endl<<copyinput.get<Matrix<double> >("mat")<<std::endl;
		std::cout<<"nthreads "<<copyinput.get<unsigned int>("nthreads")<<std::endl;
		
		//output.reset("k",278);
		//for(unsigned int i(0);i<output.size();i++){
			//std::cout<<output.name(i)<<" ";
		//}
		//std::cout<<std::endl;
		//std::cout<<output<<std::endl;
		//std::cout<<"m "<<output.get<unsigned int>("m")<<std::endl;

		GenericData<double> a("a",23.3);
		GenericData<double>* b(a.clone());
		std::cout<<b->get_val()<<std::endl;


	}
}
