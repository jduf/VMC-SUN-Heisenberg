#include "VMCMinimization.hpp"

VMCMinimization::VMCMinimization(Minimization& m, std::string const& prefix):
	m_(m),
	basename_(prefix)
{
	basename_ += "-wf" + m_.get_wf();
	basename_ += "-N"  + my::tostring(m_.get_N());
	basename_ += "-m"  + my::tostring(m_.get_m());
	basename_ += "-n"  + my::tostring(m_.get_n());
	basename_ += "-bc" + my::tostring(m_.get_bc());
}

void VMCMinimization::complete_analysis(double const& convergence_criterion){
	m_.complete_analysis(convergence_criterion);
}

void VMCMinimization::save() const { m_.save(get_filename()); }

void VMCMinimization::refine(unsigned int const& Nrefine, double const& convergence_criterion, unsigned int const& tmax){
	m_.refine(Nrefine,convergence_criterion,tmax);
}

std::shared_ptr<MCSim> VMCMinimization::compute_vmc(Vector<double> const& param){
	return m_.compute_vmc(param);
}
