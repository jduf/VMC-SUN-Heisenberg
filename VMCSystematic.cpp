#include "VMCSystematic.hpp"

VMCSystematic::VMCSystematic(VMCMinimization const& m, Vector<double> const& param, Matrix<int> const& sym, unsigned int const& p1, unsigned int const& p2):
	VMCMinimization(m,"SYS"),
	param_(param),
	sym_(sym),
	p1_(p1),
	p2_(p2)
{
	if(m_->dof_ - sym.row() != 2){
		std::cout<<__PRETTY_FUNCTION__<<"bug"<<std::endl;
	}
}

/*{public methods*/
void VMCSystematic::run(int const& nobs, double const& dE, unsigned int const& maxiter){
	for(unsigned int i(0);i<m_->ps_[p1_].size();i++){
		param_(p1_) = m_->ps_[p1_](i);
		for(unsigned int j(0);j<m_->ps_[p2_].size();j++){
			param_(p2_) = m_->ps_[p2_](j);
			apply_symmetry();
			evaluate_until_precision(param_,nobs,dE,maxiter);
		}
	}
}

void VMCSystematic::plot(){
	std::cout<<"will save "<<m_->samples_list_.size()<<" samples"<<std::endl;
	IOFiles data("systematic.dat",true);
	m_->samples_list_.set_target();
	Vector<double> param;
	while(m_->samples_list_.target_next()){
		param = m_->samples_list_.get().get_param(); 
		data<<m_->samples_list_.get().get_param()<<" "<<m_->samples_list_.get().get_MCS()->get_energy()<<IOFiles::endl;
	}

	Gnuplot gp("./","systematic");
	gp+="set xyplane 0";
	gp+="splot 'systematic.dat' u 1:9:($12+$13) w p pt 10,\\";
	gp+="      'systematic.dat' u 1:9:($12-$13) w p pt 8";
	gp.save_file();
	gp.create_image(true,true);
}

void VMCSystematic::rerun(unsigned int const& nmin, int const& nobs, double const& dE, unsigned int const& maxiter){
	double E_range(0.0);
	double tmp;
	m_->samples_list_.set_target();
	while(m_->samples_list_.target_next()){
		tmp = m_->samples_list_.get().get_MCS()->get_energy().get_x();
		if(tmp<E_range){ E_range=tmp; }
	}
	E_range *= 0.99;

	List<MCSim> sorted_list;
	m_->samples_list_.set_target();
	while(m_->samples_list_.target_next()){
		if(m_->samples_list_.get().get_MCS()->get_energy().get_x()<E_range){
			sorted_list.add_sort(m_->samples_list_.get_ptr(),MCSim::sort_by_E);
		}
	}

	sorted_list.set_target();
	unsigned int i(0);
	while(sorted_list.target_next() && i++<nmin){
		evaluate_until_precision(sorted_list.get().get_param(),nobs,dE,maxiter);
	}
}

void VMCSystematic::apply_symmetry(){
	for(unsigned int i(0);i<sym_.row();i++){
		if(sym_(i,1)<0){ param_(sym_(i,0)) = sym_(i,2); }
		else { param_(sym_(i,0)) = sym_(i,2)*param_(sym_(i,1)); }
	}
}
/*}*/
