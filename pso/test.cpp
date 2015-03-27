#include"FuncPSO.hpp"
#include"Parseur.hpp"
#include"Time.hpp"

int main(int argc,char* argv[]){
	Parseur P(argc,argv);
	Time t;
	FuncPSO s(8,100,2,2.1,2.1);

	s.init(200.0);
	//std::cout<<"initialization done in "<<t.elapsed()<<std::endl;
	//t.set();
	//s.run();
	//std::cout<<"run done in "<<t.elapsed()<<std::endl;
	//s.result();
}
