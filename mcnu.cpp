/*!  @file mcnu.cpp */

#include "Parseur.hpp"
#include "MonteCarlo.hpp"
#include "SquareJastrow.hpp"

void copy_input_mcnu(Container const& c, Container& new_c);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	//FileParser file(P.get<std::string>("nu"));
	if(!P.status()){
		Container prop;
		Container input;
		Matrix<double> nu(4,4);
		//for(unsigned int i(0);i<4;i++){
			//nu(i,0) = -0;
			//nu(i,1) = -0.269;
			//nu(i,2) = -0.269;
			//nu(i,3) = -0;
		//}
		for(unsigned int i(0);i<4;i++){
			nu(i,0) = 1.10845 ;
			nu(i,1) = 1.29585;
			nu(i,2) =-1.1572;
			nu(i,3) =-0.434311;
		}
		input.set("nu",nu);
		//file.extract<Matrix<double> >("nu",input);

		if(P.check("t_max")){ prop.set("t_max",P.get<unsigned int>("t_max")); } 
		else { unsigned int t(300); prop.set("t_max",t);}

		SquareJastrow s(P);
		s.properties(prop);
		copy_input_mcnu(prop,input);

		MonteCarlo<std::complex<double> > sim(input);
		Matrix<unsigned int> lattice;
		sim.run(3,1e6,&lattice);
		s.lattice(lattice);
		std::cout<<lattice<<std::endl;
		std::cout<<"energy = "<<sim.get_energy()<<std::endl;;
	}
}

void copy_input_mcnu(Container const& c, Container& new_c){
	new_c.set("fermionic",false);
	new_c.set("t_max",c.get<unsigned int>("t_max"));
	new_c.set("N",c.get<unsigned int>("N"));
	new_c.set("m",c.get<unsigned int>("m"));
	new_c.set("sl",c.get<Vector<unsigned int> >("sl"));
	new_c.set("nn",c.get<Matrix<unsigned int> >("nn"));
	new_c.set("cc",c.get<Matrix<unsigned int> >("cc"));
	new_c.set("sts",c.get<Matrix<unsigned int> >("sts"));
	new_c.set("omega",c.get<Matrix<std::complex<double> > >("omega"));
}
