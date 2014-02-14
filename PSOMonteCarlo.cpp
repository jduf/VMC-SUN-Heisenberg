#include "PSOMonteCarlo.hpp"

PSOMonteCarlo::PSOMonteCarlo(Parseur& P,unsigned int Nfreedom):
	PSO(P.get<unsigned int>("Nbees"),Nfreedom,P.get<double>("cg"),P.get<double>("cp"),P.get<double>("maxiter")),
	input_(false),
	param_(true),
	results_("bla"),
	z_(0)
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

PSOMonteCarlo::~PSOMonteCarlo(){
	Write nu("nu.jdbin");
	nu<<create_nu(pb_[bbee_]);
}

void PSOMonteCarlo::get_system_properties(Parseur& P, Container& c){
	wf_=P.get<std::string>("wf");
	if( wf_ == "trianglejastrow" ){
		TriangleJastrow s(P);
		s.properties(c);
		z_=6;
	}
	if( wf_ == "jastrow" ){
		SquareJastrow s(P);
		s.properties(c);
		z_=4;
	}
}

void PSOMonteCarlo::copy_input(Container const& c, Container& new_c){
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

double PSOMonteCarlo::run(Vector<double> const& x){
	Container input_new;
	copy_input(input_,input_new);

	input_new.set("nu",create_nu(x));
	MonteCarlo<std::complex<double> > sim(input_new);
	sim.run(2,1e5);
	std::cout<<x<<" "<<sim.get_energy()<<" "<<sim.get_status()<<std::endl;
	return sim.get_energy();
}

Matrix<double> PSOMonteCarlo::create_nu(Vector<double> const& x){
	Matrix<double> nu(z_,Nfreedom_+1);
	for(unsigned int i(0);i<z_;i++){
		for(unsigned int j(0);j<Nfreedom_;j++){
			nu(i,j) = x(j);
		}
		nu(i,Nfreedom_) = 0;
	}
	return nu;
}

void PSOMonteCarlo::save(){
	for(unsigned int i(0);i<Nbees_;i++){
		results_<<pb_[i]<<" "<<pfb_[i]<<Write::endl;
	}
}
