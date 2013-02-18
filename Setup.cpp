#include "Write.hpp"
#include "Read.hpp"
#include "Parseur.hpp"
#include <sstream>

//int main(int argc, char* argv[]){
	//Parseur P(argc,argv);
	//unsigned int N_spin(4);
	//unsigned int N_m(6);
	//unsigned int N_n(0);
	//std::string matrice;
	//std::string lattice;
	//std::string sysname;
//
	//P.set("N_spin",N_spin);	
	//P.set("N_m",N_m);	
	//P.set("N_n",N_n);	
	//P.set("matrice",matrice);	
//
	//std::stringstream ss1;
	//std::stringstream ss2;
	//ss1<<N_spin;
	//ss2<<N_spin*N_m;
	//switch(N_n){
		//case 2:
			//sysname = "chain";
			//break;
		//case 3:
			//sysname = "honeycomb";
			//break;
		//case 4:
			//sysname =  "square";
			//break;
	//}
	//sysname += "-N"+ss1.str() + "-S"+ss2.str()+".bin";
//
	//Matrice<double> m(N_spin*N_m);
	//Read r(matrice.c_str(),false);
	//r>>m;
	//Write w(sysname.c_str());
	//w<<N_spin<<N_m<<N_n<<m;
//
	//unsigned int N_spin_out(0);
	//unsigned int N_m_out(0);
	//unsigned int N_n_out(0);
	//Read r_out(sysname.c_str());
	//Matrice<double> U;
	//r_out>>N_spin_out>>N_m_out>>N_n_out>>U;
	//std::cout<<N_spin_out<<" "<<N_m_out<<" "<<N_n_out<<std::endl;
	//U.print();
//}
int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::string sysname;
	unsigned int N_spin(0), N_m(0), N_n(0);
	P.set("sysname",sysname);	

	Read r(sysname.c_str());
	r>>N_spin>>N_m>>N_n;
	Matrice<double> m(N_spin*N_m);
	r>>m;

	m.print();
	Write w(sysname.c_str());
	bool need_complex(false);
	w<<N_spin<<N_m<<N_n<<need_complex<<m;
}
