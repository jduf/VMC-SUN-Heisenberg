#include "FuncPSO.hpp"

FuncPSO::FuncPSO(unsigned int Nparticles, unsigned int maxiter, unsigned int Nfreedom, double cg, double cp):
	Swarm<ParticleOnGrid>(Nparticles,maxiter,Nfreedom,cg,cp)
{
	for(unsigned int i(0);i<Nfreedom;i++){
		Particle::set_limit(i,0,2);
	}
}

double FuncPSO::f(Vector<double> const& x){
	return (1-x(0)*x(0))*(1-x(0)*x(0))+(3-x(1)*x(1))*(3-x(1)*x(1));
}

bool FuncPSO::evaluate(unsigned int const& p){
	std::shared_ptr<ParticleOnGrid> particle(std::dynamic_pointer_cast<ParticleOnGrid>(particle_[p]));
	double fx(f(particle->get_x()));
	std::shared_ptr<Measure> result(std::make_shared<Measure>(particle->get_x(),fx));
#pragma omp critical
	{
		m_.add_sort(result,Measure::func);
	}

	if( fx<particle_[p]->get_fbx() ){ 
		particle->update(fx);
		return true;
	} else {
		return false;
	}
}

void FuncPSO::result(){
	std::cout<<m_<<std::endl;
	std::cout<<"best for each particle"<<std::endl;
	print();
}

