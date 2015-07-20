#include"FuncPSO.hpp"
//#include"Parseur.hpp"
//#include"Time.hpp"

//int main(int argc,char* argv[]){
	//Parseur P(argc,argv);
int main(){
	FuncPSO s(16,100,2,2.1,2.1);

	s.init_PSO(200.0);
	s.minimize();
}
