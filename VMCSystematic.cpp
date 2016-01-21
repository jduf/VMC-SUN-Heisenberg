#include "VMCSystematic.hpp"

VMCSystematic::VMCSystematic(VMCMinimization const& m, Parseur& P):
	VMCMinimization(m,"SYS")
{
	set_phase_space(P);
}

/*{public methods*/
void VMCSystematic::run(int const& nobs, double const& dE, unsigned int const& maxiter){
	Vector<double> param(1);
	for(unsigned int i(0);i<m_->ps_[0].size();i++){
		param(0) = m_->ps_[0](i);
		std::cout<<std::endl;
		std::cout<<param<<std::endl;
		evaluate_until_precision(param,nobs,dE,maxiter);
	}
}

void VMCSystematic::plot(){
	std::cout<<"will save "<<m_->samples_.size()<<" samples"<<std::endl;
	IOFiles data("systematic.dat",true);
	m_->samples_.set_target();
	Vector<double> param;
	while(m_->samples_.target_next()){
		param = m_->samples_.get().get_param(); 
		data<<m_->samples_.get().get_param()<<" "<<m_->samples_.get().get_MCS()->get_energy()<<IOFiles::endl;
	}

	Gnuplot gp("./","systematic");
	gp+="set xyplane 0";
	gp+="plot 'systematic.dat' u 1:2:3 w e notitle";
	gp.save_file();
	gp.create_image(true,true);
}
/*}*/
