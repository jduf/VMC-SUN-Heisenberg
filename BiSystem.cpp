#include "BiSystem.hpp"

/*constructors, destructors*/
/*{*/
BiSystem::BiSystem(Parseur& P):
	s_(P)
{}

BiSystem::BiSystem(IOFiles& r):
	s_(r)
{
	unsigned int size(r.read<unsigned int>());
	for(unsigned int i(0);i<size;i++){ param_.push_back(Vector<double>(r)); }
	for(unsigned int i(0);i<size;i++){
		mcsys_.push_back(std::vector<std::unique_ptr<MCSystem> >());
		for(unsigned int j(0);j<size;j++){
			if(i!=j){ mcsys_[i].push_back(std::unique_ptr<MCSystem>(new SystemBiFermionic<double>(r))); }
			else { mcsys_[i].push_back(std::unique_ptr<MCSystem>(new SystemFermionic<double>(r))); }
		}
	}
	r>>H_>>O_>>dH_>>dO_>>E_>>dE_;
}
/*}*/

/*core methods*/
/*{*/
void BiSystem::add_new_param(Vector<double> const& param){
	std::cout<<param<<std::endl;
	CreateSystem cs(&s_);
	cs.init(&param,NULL);
	if(cs.get_status()==2){
		cs.create(true);
		if(cs.get_status()==1){
			param_.push_back(param);

			mcsys_.push_back(std::vector<std::unique_ptr<MCSystem> >(param_.size()));
			for(unsigned int i(0);i<param_.size();i++){ mcsys_[i].push_back(std::unique_ptr<MCSystem>() ); }

			unsigned int j(param_.size()-1);
			if(cs.use_complex()){
				mcsys_[j][j].reset(new SystemFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(cs.get_GenericSystem())));
			} else {
				mcsys_[j][j].reset(new SystemFermionic<double>(*dynamic_cast<const Fermionic<double>*>(cs.get_GenericSystem())));
			}

			for(unsigned int i(0);i<j;i++){
				if(cs.use_complex()){
					mcsys_[i][j].reset(new SystemBiFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(mcsys_[i][i].get()),*dynamic_cast<const Fermionic<std::complex<double> >*>(mcsys_[j][j].get())));
					mcsys_[j][i].reset(new SystemBiFermionic<std::complex<double> >(*dynamic_cast<const Fermionic<std::complex<double> >*>(mcsys_[j][j].get()),*dynamic_cast<const Fermionic<std::complex<double> >*>(mcsys_[i][i].get())));
				} else {
					mcsys_[i][j].reset(new SystemBiFermionic<double>(*dynamic_cast<const Fermionic<double>*>(mcsys_[i][i].get()),*dynamic_cast<const Fermionic<double>*>(mcsys_[j][j].get())));
					mcsys_[j][i].reset(new SystemBiFermionic<double>(*dynamic_cast<const Fermionic<double>*>(mcsys_[j][j].get()),*dynamic_cast<const Fermionic<double>*>(mcsys_[i][i].get())));
				}
			}
		} else { std::cout<<__PRETTY_FUNCTION__<<" : CreateSystem::create(&p,NULL) failed "<<std::endl; }
	} else { std::cout<<__PRETTY_FUNCTION__<<" : CreateSystem::init(&p,NULL) failed "<<std::endl; }
}

void BiSystem::run(unsigned int const& nruns, unsigned int const& tmax){
	for(unsigned int i(0);i<param_.size();i++){
		for(unsigned int j(0);j<param_.size();j++){
#pragma omp parallel for
			for(unsigned int k=0;k<nruns;k++){
				std::unique_ptr<MCSystem> tmp(mcsys_[i][j]->clone());

				MonteCarlo sim(tmp.get(),(i==j?tmax:50*tmax));
				sim.thermalize(1e6);
				sim.run();

#pragma omp critical(System__merge)
				{ mcsys_[i][j]->merge(tmp.get()); }
			}
			mcsys_[i][j]->complete_analysis(1e-5);
			mcsys_[i][j]->print(2);
		}
	}
}

void BiSystem::compute_E(){
	unsigned int n_wfs(param_.size());// # of wavefunctions
	H_.set(n_wfs,n_wfs);
	dH_.set(n_wfs,n_wfs);
	O_.set(n_wfs,n_wfs);
	dO_.set(n_wfs,n_wfs);
	E_.set(n_wfs);
	dE_.set(n_wfs);

	for(unsigned int i(0);i<n_wfs;i++){
		for(unsigned int j(0);j<n_wfs;j++){
			H_(i,j) = mcsys_[i][j]->get_energy().get_x();
			dH_(i,j)= mcsys_[i][j]->get_energy().get_dx();
			if(i==j){
				O_(i,j) = 1;
				dO_(i,j)= 0;
			} else {
				O_(i,j) = mcsys_[i][j]->get_obs()[1][0].get_x();
				dO_(i,j)= mcsys_[i][j]->get_obs()[1][0].get_dx();
			}
		}
	}

	for(unsigned int i(0);i<n_wfs;i++){
		for(unsigned int j(i+1);j<n_wfs;j++){
			/*be careful with complex values, not to scew the phase taking the square root*/
			H_(i,j) = sqrt(H_(i,j)*H_(j,i));
			H_(i,j) = H_(j,i);
			dH_(i,j)= (dH_(i,j)+dH_(j,i))/(2*sqrt(2));
			dH_(i,j)= dH_(j,i);
			O_(i,j) = sqrt(O_(i,j)*O_(j,i));
			O_(i,j) = O_(j,i);
			dO_(i,j)= (dO_(i,j)+dO_(j,i))/(2*sqrt(2));
			dO_(i,j)= dO_(j,i);
		}
	}

	Lapack<double>(H_,true,'S').generalized_eigensystem(O_,E_);
}

void BiSystem::compute_dE(){
	unsigned int n_wfs(param_.size());// # of wavefunctions
	unsigned int n_rnd(1000); 	 	  // # of random evaluations
	unsigned int n_bin(20);		 	  // # of bin for the history

	Matrix<double> Htmp;
	Matrix<double> Otmp;
	Vector<double> Etmp;

	Matrix<double> EM(n_wfs,n_rnd); // the eigenvalues are stored in vertical

	Rand<double> rnd(-1.0,1.0);
	for(unsigned int i(0);i<n_rnd;i++){
		Htmp = H_;
		Otmp = O_;
		Etmp.set();
		for(unsigned int j(0);j<n_wfs;j++){
			for(unsigned int k(j);k<n_wfs;k++){
				Htmp(j,k)+= rnd.get()*dH_(j,k);
				Htmp(k,j) = Htmp(j,k);
				Otmp(j,k)+= rnd.get()*dO_(j,k);
				Otmp(k,j) = Otmp(j,k);
			}
		}
		Lapack<double>(Htmp,false,'S').generalized_eigensystem(Otmp,Etmp);
		if(Etmp.ptr()){ for(unsigned int j(0);j<n_wfs;j++){ EM(j,i) = Etmp(j); } }
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
		dE_(i) = p(1);
	}
}
/*}*/

void BiSystem::save() const {
	CreateSystem cs(&s_);
	cs.init(&param_[0],NULL);

	IOFiles w(cs.get_filename()+".jdbin",true);
	s_.write(w);
	unsigned int size(param_.size());
	w<<size;
	for(unsigned int i(0);i<size;i++){ w<<param_[i]; }
	for(unsigned int i(0);i<size;i++){
		for(unsigned int j(0);j<size;j++){
			mcsys_[i][j]->write(w);
		}
	}
	w<<H_<<O_<<dH_<<dO_<<E_<<dE_;
}
