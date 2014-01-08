#ifndef DEF_PSOMONTECARLO
#define DEF_PSOMONTECARLO

#include "MonteCarlo.hpp"
#include "Parseur.hpp"
#include "SquareJastrow.hpp"
#include "TriangleJastrow.hpp"
#include "PSO.hpp"

class PSOMonteCarlo : public PSO {
	public:
		PSOMonteCarlo(Parseur& P);

		void save();

	private:
		virtual double run(Vector<double> const& x);
		void copy_input(Container const& input, Container& new_input);
		void get_system_properties(Parseur& P, Container& c);

		Container input_;
		Container param_;
		Write results_;
		unsigned int iter;
		std::string wf_;
};

PSOMonteCarlo::PSOMonteCarlo(Parseur& P):
	PSO(P.get<unsigned int>("Nbees"),P.get<unsigned int>("Nfreedom"),P.get<double>("cg"),P.get<double>("cp"),P.get<double>("maxiter")),
	input_(false),
	param_(true),
	results_("bla"),
	iter(0)
{
	if(!P.status()){
		Container prop;
		get_system_properties(P,prop);

		if(P.check("t_max")){ prop.set("t_max",P.get<unsigned int>("t_max")); } 
		else { unsigned int t(300); prop.set("t_max",t);}

		copy_input(prop,input_);

		param_.set("Lx",prop.get<unsigned int>("Lx"));
		param_.set("Ly",prop.get<unsigned int>("Ly"));
		param_.set("bc",prop.get<double>("bc"));

	} else {
		std::cerr<<"blof"<<std::endl;
	}
}

void PSOMonteCarlo::get_system_properties(Parseur& P, Container& c){
	wf_=P.get<std::string>("wf");
	if( wf_ == "trianglejastrow" ){
		TriangleJastrow s(P);
		s.properties(c);
	}
	if( wf_ == "jastrow" ){
		SquareJastrow s(P);
		s.properties(c);
	}
}

void PSOMonteCarlo::copy_input(Container const& c, Container& new_c){
	new_c.set("fermionic",false);
	new_c.set("t_max",c.get<unsigned int>("t_max"));
	new_c.set("N",c.get<unsigned int>("N"));
	new_c.set("m",c.get<unsigned int>("m"));
	new_c.set("sl",c.get<Vector<unsigned int> >("sl"));
	new_c.set("nn",c.get<Matrix<unsigned int> >("nn"));
	new_c.set("sts",c.get<Matrix<unsigned int> >("sts"));
	new_c.set("omega",c.get<Matrix<std::complex<double> > >("omega"));
}

double PSOMonteCarlo::run(Vector<double> const& x){
	Container input_new;
	copy_input(input_,input_new);
	Vector<double> nu;
	if( wf_ == "jastrow" ){
		nu.set(12);
		nu(0)=x(0);
		nu(1)=x(1);
		nu(2)=x(0);
		nu(3)=x(1);
		nu(4)=x(2);
		nu(5)=x(3);
		nu(6)=x(4);
		nu(7)=x(5);
		nu(8)=x(2);
		nu(9)=x(3);
		nu(10)=x(4);
		nu(11)=x(5);
	}
	if( wf_ == "trianglejastrow" ){
		nu.set(18);
		nu(0)=x(0);
		nu(1)=x(1);
		nu(2)=x(2);
		nu(3)=x(0);
		nu(4)=x(1);
		nu(5)=x(2);
		nu(6)=x(3);
		nu(7)=x(4);
		nu(8)=x(5);
		nu(9)=x(6);
		nu(10)=x(7);
		nu(11)=x(8);
		nu(12)=x(3);
		nu(13)=x(4);
		nu(14)=x(5);
		nu(15)=x(6);
		nu(16)=x(7);
		nu(17)=x(8);
	}
	
	input_new.set("nu",nu);
	input_new.set("filename","data"+tostring(iter++));
	MonteCarlo<std::complex<double> > sim(input_new);
	sim.run(2,1e6);
	std::cout<<x<<" "<<sim.get_energy()<<std::endl;
	return sim.get_energy();
}

void PSOMonteCarlo::save(){
	results_<<"%";
	for(unsigned int i(0);i<Nbees_;i++){
		results_<<pb_[i]<<" "<<pfb_[i]<<Write::endl;
	}
}
#endif

