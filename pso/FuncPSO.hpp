#ifndef DEF_FUNCPSO
#define DEF_FUNCPSO

#include "PSO.hpp"
#include "List.hpp"

class Measure{
	public:
		Measure(Vector<double> const& x, double const& fx):fx_(fx),x_(x),N_(0){};

		void print(std::ostream& flux) const
		{ flux<<x_<<" : "<<fx_<<std::endl; }

		double fx_;
		Vector<double> x_;
		unsigned int N_;
};

bool func(Measure const& a, Measure const& b) { 
	unsigned int i(0);
	while(i<a.x_.size()){
		if(a.x_(i) > b.x_(i)){ return false; }
		if(a.x_(i) < b.x_(i)){ return true; }
		if(a.x_(i)== b.x_(i)){ i++; }
	}
	return false;
}

class ParticleOnGrid: public Particle {
	public:
		ParticleOnGrid(){};
		~ParticleOnGrid(){
#pragma omp critical
			{
				std::cout<<particle_history_<<std::endl;
			}
		}
		void move(Vector<double> const& bx_all);

		void add_history(std::shared_ptr<Measure> const& m);

	protected:
		List<Measure> particle_history_;
};

void ParticleOnGrid::move(Vector<double> const& bx_all){
	Particle::move(bx_all);
	unsigned int n;
	double dx(0.01);
	for(unsigned int j(0);j<Nfreedom_;j++){
		n=0;
		if(std::abs(x_(j))<dx/2){ n=1; }
		if(std::abs(x_(j)-min_(j))<dx/2){ n=2; }
		if(std::abs(x_(j)-max_(j))<dx/2){ n=3; }
		switch(n){
			case 0:{ x_(j) = std::round(x_(j)/dx)*dx; }break;
			case 1:{ x_(j) = 0; }break;
			case 2:{ x_(j) = min_(j); }break;
			case 3:{ x_(j) = max_(j); }break;
		}
	}
}

void ParticleOnGrid::add_history(std::shared_ptr<Measure> const& m){
	particle_history_.add_sort(m,func);
}

std::ostream& operator<<(std::ostream& flux, Measure const& m){
	m.print(flux);
	return flux;
}

class FuncPSO : public Swarm<ParticleOnGrid> {
	public:
		FuncPSO(unsigned int Nparticles, unsigned int maxiter, unsigned int Nfreedom, double cg, double cp);
		~FuncPSO(){};
		void result();

		double f(Vector<double> const& x);

	protected:
		List<Measure> m_;

	private:
		bool evaluate(unsigned int const& p); 
};

FuncPSO::FuncPSO(unsigned int Nparticles, unsigned int maxiter, unsigned int Nfreedom, double cg, double cp):
	Swarm<ParticleOnGrid>(Nparticles,maxiter,Nfreedom,cg,cp)
{}

double FuncPSO::f(Vector<double> const& x){
	return (1-x(0)*x(0))*(1-x(0)*x(0))+(3-x(1)*x(1))*(3-x(1)*x(1));
}

bool FuncPSO::evaluate(unsigned int const& p){
	double fx(f(particle_[p]->get_x()));
	std::shared_ptr<Measure> result(std::make_shared<Measure>(particle_[p]->get_x(),fx));
#pragma omp critical
	{
		m_.add_sort(result,func);
	}

	if( fx<particle_[p]->get_fbx() ){ 
		particle_[p]->set_best(particle_[p]->get_x(),fx);
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
#endif
