/*!@file mcbi.cpp */

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"
#include <omp.h>

void run(Matrix<double>& H, Matrix<double>& dH, Matrix<double>& O, Matrix<double>& dO, Matrix<double> const& MP, System const& sys, unsigned int const& tmax, unsigned int const& nruns, unsigned int const& i, unsigned int const& j);
void run(Matrix<double>& H, Matrix<double>& dH, Matrix<double>& O, Matrix<double>& dO, Matrix<double> const& MP, System const& sys, unsigned int const& tmax, unsigned int const& nruns, unsigned int const& i);
void compute_E(Matrix<double>& H, Matrix<double>& dH, Matrix<double>& O, Matrix<double>& dO);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int i(0);
	unsigned int tmax(P.get<unsigned int>("tmax"));
	unsigned int nruns(P.find("nruns",i,false)?P.get<unsigned int>(i):omp_get_max_threads());
	if(!P.find("M",i,false)){
		std::vector<unsigned int> M(P.get<unsigned int>("N"),P.get<unsigned int>("n")*P.get<unsigned int>("m")/P.get<unsigned int>("N"));
		P.set("M",M);
	}

	System sys(P);
	IOFiles r(P.get<std::string>("params"),false);
	Matrix<double> MP(P.get<unsigned int>("nwfs"),20);
	if(!P.locked()){
		r>>MP;
		std::cout<<MP<<std::endl;
		Matrix<double> H(MP.row(),MP.row(),0);
		Matrix<double> O(H);
		Matrix<double> dH(H);
		Matrix<double> dO(H);

		for(unsigned int i(0);i<MP.row();i++){
			for(unsigned int j(0);j<MP.row();j++){
				if(i!=j){ run(H,dH,O,dO,MP,sys,tmax*50,nruns,i,j); }
				else { run(H,dH,O,dO,MP,sys,tmax,nruns,i); }
			}
		}

		compute_E(H,dH,O,dO);
	} else { std::cout<<__PRETTY_FUNCTION__<<" : Parseur locked"<<std::endl; }
}

void run(Matrix<double>& H, Matrix<double>& dH, Matrix<double>& O, Matrix<double>& dO, Matrix<double> const& MP, System const& sys, unsigned int const& tmax, unsigned int const& nruns, unsigned int const& i, unsigned int const& j){
	Vector<double> p0(MP.col());
	Vector<double> p1(MP.col());
	for(unsigned int k(0);k<MP.col();k++){
		p0(k) = MP(i,k);
		p1(k) = MP(j,k);
	}

	CreateSystem cs0(&sys);
	CreateSystem cs1(&sys);
	cs0.init(&p0,NULL);
	cs1.init(&p1,NULL);
	cs0.set_obs(0);
	cs1.set_obs(0);
	if(cs0.get_status()==2 && cs1.get_status()==2){
		cs0.create(true);
		cs1.create(true);
		if(cs0.get_status()==1 && cs1.get_status()==1){
#pragma omp parallel for
			for(unsigned int j=0;j<nruns;j++){
				MCSystem* mcsys(NULL);
				if(cs0.use_complex() && cs1.use_complex()){
					mcsys = new SystemBiFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs0.get_GenericSystem()),*dynamic_cast<const Fermionic<std::complex<double> >*>(cs1.get_GenericSystem()));
				}
				if(!cs0.use_complex() && !cs1.use_complex()){
					mcsys = new SystemBiFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs0.get_GenericSystem()),*dynamic_cast<const Fermionic<double>*>(cs1.get_GenericSystem()));
				}
				if(!mcsys){ std::cout<<__PRETTY_FUNCTION__<<" MCSystem was not constructed"<<std::endl; }
				else {
					MonteCarlo sim(mcsys,tmax);
					sim.thermalize(1e6);
					sim.run(1e7);

#pragma omp critical(System__merge)
					{ cs0.merge(mcsys); }

					delete mcsys;
				}
			}
			cs0.complete_analysis(1e-5);
			H(i,j) = cs0.get_obs()[0][0].get_x();
			dH(i,j)= cs0.get_obs()[0][0].get_dx();
			O(i,j) = cs0.get_obs().back()[0].get_x();
			dO(i,j)= cs0.get_obs().back()[0].get_dx();
		} else { std::cout<<__PRETTY_FUNCTION__<<" : CreateSystem::create(&p,NULL) failed "<<std::endl; }
	} else { std::cout<<__PRETTY_FUNCTION__<<" : CreateSystem::init(&p,NULL) failed "<<std::endl; }
}

void run(Matrix<double>& H, Matrix<double>& dH, Matrix<double>& O, Matrix<double>& dO, Matrix<double> const& MP, System const& sys, unsigned int const& tmax, unsigned int const& nruns, unsigned int const& i){
	Vector<double> p(MP.col());
	for(unsigned int k(0);k<MP.col();k++){
		p(k) = MP(i,k);
	}

	CreateSystem cs(&sys);
	cs.init(&p,NULL);
	cs.set_obs(0);
	if(cs.get_status()==2){
		cs.create(true);
		if(cs.get_status()==1){
#pragma omp parallel for
			for(unsigned int j=0;j<nruns;j++){
				MCSystem* mcsys(NULL);
				if(cs.use_complex()){
					mcsys = new SystemFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs.get_GenericSystem()));
				}
				if(!cs.use_complex()){
					mcsys = new SystemFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs.get_GenericSystem()));
				}
				if(!mcsys){ std::cout<<__PRETTY_FUNCTION__<<" MCSystem was not constructed"<<std::endl; }
				else {
					MonteCarlo sim(mcsys,tmax);
					sim.thermalize(1e6);
					sim.run(1e8);

#pragma omp critical(System__merge)
					{ cs.merge(mcsys); }

					delete mcsys;
				}
			}
			cs.complete_analysis(1e-5);
			H(i,i) = cs.get_obs()[0][0].get_x();
			dH(i,i) = cs.get_obs()[0][0].get_dx();
			O(i,i) = 1;
			dO(i,i) = 0;
		} else { std::cout<<__PRETTY_FUNCTION__<<" : CreateSystem::create(&p,NULL) failed "<<std::endl; }
	} else { std::cout<<__PRETTY_FUNCTION__<<" : CreateSystem::init(&p,NULL) failed "<<std::endl; }
}

void compute_E(Matrix<double>& H, Matrix<double>& dH, Matrix<double>& O, Matrix<double>& dO){
	unsigned int n_wfs(H.row());// number of wavefunctions
	unsigned int n_rnd(1000);	// number of random evaluations
	unsigned int n_bin(20);		// number of bin for the history

	for(unsigned int i(0);i<n_wfs;i++){
		for(unsigned int j(i+1);j<n_wfs;j++){
			H(i,j) = sqrt(H(i,j)*H(j,i));
			H(i,j) = H(j,i);
			dH(i,j)= (dH(i,j)+dH(j,i))/(2*sqrt(2));
			dH(i,j)= dH(j,i);
			O(i,j) = sqrt(O(i,j)*O(j,i));
			O(i,j) = O(j,i);
			dO(i,j)= (dO(i,j)+dO(j,i))/(2*sqrt(2));
			dO(i,j)= dO(j,i);
		}
	}

	Matrix<double> Htmp;
	Matrix<double> Otmp;

	Vector<double> E(n_wfs);
	Vector<double> dE(n_wfs);
	Matrix<double> EM(n_wfs,n_rnd); // the eigenvalues are stored in vertical

	Rand<double> rnd(-1,1);
	for(unsigned int i(0);i<n_rnd;i++){
		Htmp = H;
		Otmp = O;
		for(unsigned int j(0);j<n_wfs;j++){
			for(unsigned int k(j);k<n_wfs;k++){
				Htmp(j,k)+= rnd.get()*dH(j,k);
				Htmp(k,j) = Htmp(j,k);
				Otmp(j,k)+= rnd.get()*dO(j,k);
				Otmp(k,j) = Otmp(j,k);
			}
		}
		Lapack<double>(Htmp,false,'S').generalized_eigensystem(Otmp,E);
		for(unsigned int j(0);j<n_wfs;j++){ EM(j,i) = E(j); }
	}

	double min;
	double max;
	double dx;

	Vector<double> x(n_bin);
	Vector<double> y(n_bin,0.0);

	auto gaussian = [](double const& x, const double* p){
		return exp(-(x-p[0])*(x-p[0])/(2*p[1]))/p[2];
	};

	for(unsigned int i(0);i<n_wfs;i++){
		min = EM(i,0);
		max = EM(i,0);
		for(unsigned int j(1);j<n_rnd;j++){
			if( EM(i,j) < min ){ min = EM(i,j); }
			if( EM(i,j) > max ){ max = EM(i,j); }
		}
		dx = (max-min)/n_bin;
		for(unsigned int j(0);j<n_bin;j++){ 
			x(j) = min+j*dx; 
			y(j) = 0.0;
		}
		for(unsigned int j(0);j<n_rnd;j++){
			for(unsigned int k(0);k<n_bin;k++){
				if(EM(i,j)>=x(k)){ y(k)+=1.0; k=n_bin; }
			}
		}
		y /= n_rnd;

		Vector<double> p(3);
		p(0) = x.mean();
		p(1) = x.variance();
		p(2) = 1;
		Fit(x,y,p,gaussian);
		dE(i) = p(1);
	}

	Lapack<double>(H,true,'S').generalized_eigensystem(O,E);

	//std::cout<<H<<std::endl;
	//std::cout<<std::endl;
	//std::cout<<O<<std::endl;
	//H.print_mathematica();
	//O.print_mathematica();
	
	std::cout<<E<<std::endl;
	std::cout<<dE<<std::endl;
}
