/*!  @file mcbi.cpp */

#include "MonteCarlo.hpp"
#include "CreateSystem.hpp"
#include <omp.h>

void run(Matrix<double>& H, Matrix<double>& O, Matrix<double> const& MP, System const& sys, unsigned int const& tmax, unsigned int const& nruns, unsigned int const& i, unsigned int const& j);
void run(Matrix<double>& H, Matrix<double>& O, Matrix<double> const& MP, System const& sys, unsigned int const& tmax, unsigned int const& nruns, unsigned int const& i);

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
	if(!P.locked()){
		IOFiles r("param.dat",false);
		Matrix<double> MP(2,18);
		r>>MP;
		std::cout<<MP<<std::endl;
		Matrix<double> H(MP.row(),MP.row(),0);
		Matrix<double> O(MP.row(),MP.row(),0);

		for(unsigned int i(0);i<MP.row();i++){
			for(unsigned int j(0);j<MP.row();j++){
				if(i!=j){ run(H,O,MP,sys,tmax,nruns,i,j); }
				else { run(H,O,MP,sys,tmax,nruns,i); }
			}
		}

		for(unsigned int i(0);i<MP.row();i++){
			for(unsigned int j(i+1);j<MP.row();j++){
				H(i,j) = sqrt(H(i,j)*H(j,i));
				H(i,j) = H(j,i);
				O(i,j) = sqrt(O(i,j)*O(j,i));
				O(i,j) = O(j,i);
			}
		}

		std::cout<<H<<std::endl;
		std::cout<<std::endl;
		std::cout<<O<<std::endl;

		H.print_mathematica();
		O.print_mathematica();

		Vector<double> E;
		Lapack<double>(H,true,'S').generalized_eigensystem(O,E);
		std::cout<<E<<std::endl;

	} else { std::cout<<__PRETTY_FUNCTION__<<" : Parseur locked"<<std::endl; }
}

void run(Matrix<double>& H, Matrix<double>& O, Matrix<double> const& MP, System const& sys, unsigned int const& tmax, unsigned int const& nruns, unsigned int const& i, unsigned int const& j){
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
					MonteCarlo sim(mcsys,tmax*10);
					sim.thermalize(1e6);
					sim.run(1e7);

#pragma omp critical(System__merge)
					{ cs0.merge(mcsys); }

					delete mcsys;
				}
			}
			cs0.complete_analysis(1e-5);
			H(i,j) = cs0.get_obs()[0][0].get_x();
			O(i,j) = cs0.get_obs().back()[0].get_x();
		} else { std::cout<<__PRETTY_FUNCTION__<<" : CreateSystem::create(&p,NULL) failed "<<std::endl; }
	} else { std::cout<<__PRETTY_FUNCTION__<<" : CreateSystem::init(&p,NULL) failed "<<std::endl; }
}

void run(Matrix<double>& H, Matrix<double>& O, Matrix<double> const& MP, System const& sys, unsigned int const& tmax, unsigned int const& nruns, unsigned int const& i){
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
			O(i,i) = 1;
		} else { std::cout<<__PRETTY_FUNCTION__<<" : CreateSystem::create(&p,NULL) failed "<<std::endl; }
	} else { std::cout<<__PRETTY_FUNCTION__<<" : CreateSystem::init(&p,NULL) failed "<<std::endl; }
}
