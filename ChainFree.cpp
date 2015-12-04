#include "ChainFree.hpp"

ChainFree::ChainFree(System const& s, Vector<double> const& t, Vector<double> const& mu):
	System(s),
	Chain<double>(set_spuc(t,mu,N_/m_),"chain-free"),
	t_(t),
	mu_(mu)
{
	if(status_==2){
		init_fermionic();
		filename_ += "-t";
		for(unsigned int i(0);i<t_.size();i++){ filename_ += ((t_(i)>0)?"+":"")+my::tostring(t_(i)); }
		filename_ += "-mu";
		for(unsigned int i(0);i<mu_.size();i++){ filename_ += ((mu_(i)>0)?"+":"")+my::tostring(mu_(i)); }

		system_info_.text("Trial wavefunction with different real hopping and chemical potential");
	}
}

/*{method needed for running*/
void ChainFree::compute_H(unsigned int const& c){
	(void)(c);
	H_.set(n_,n_,0);

	unsigned int s0(0);
	unsigned int s1(0);
	unsigned int ab(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		ab = s0%spuc_;
		H_(s0,s1) = obs_[0](i,4)*t_(ab);
		//if(ab==c){ H_(s0,s0) = mu_(ab); }
		H_(s0,s0) = mu_(ab);
	}
	H_ += H_.transpose();
}

void ChainFree::create(){
	for(unsigned int c(0);c<N_;c++){
		compute_H(c);
		diagonalize(true);
		if(status_==1){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		}
		if(c!=N_-1){status_++;}
	}
}

void ChainFree::save_param(IOFiles& w) const {
	std::string t_string("");
	for(unsigned int i(0);i<t_.size()-1;i++){ t_string += my::tostring(t_(i))+","; }
	t_string += my::tostring(t_.back());
	w.write("t ("+t_string+")",t_);

	std::string mu_string("");
	for(unsigned int i(0);i<mu_.size()-1;i++){ mu_string += my::tostring(mu_(i))+","; }
	mu_string += my::tostring(mu_.back());
	w.write("mu ("+mu_string+")",mu_);
}

unsigned int ChainFree::set_spuc(Vector<double> const& t, Vector<double> const& mu, unsigned int const& spuc){
	if(t.size() == spuc && mu.size() == spuc && !my::are_equal(t,Vector<double>(spuc,1.0))){ return spuc; }
	else {
		std::cerr<<__PRETTY_FUNCTION__<<" : invalid or incoherent t and mu sizes : t:="<<t.size()<<", mu:="<<mu.size()<<std::endl;
		return spuc+1;
	}
}
/*}*/

/*{method needed for checking*/
void ChainFree::check(){
	long_range_correlation_and_structure_factor();
}

void ChainFree::energy_bound(){
	compute_H(0);

	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	double x_shift(spuc_+2);

	PSTricks ps(info_+path_+dir_,filename_+"-pstricks");
	ps.begin(-1,-5,n_/1.5,2,filename_+"-pstricks");
	double t;
	double mu;
	double corr;
	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<spuc_;i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		xy0(0) = s0;
		xy1(0) = s1;

		t=H_(s0,s1);
		if(std::abs(t)>1e-4){
			if(xy1(0)<xy0(0)){
				xy1(0) = xy0(0)+1;
				linestyle="dashed";
			} else { linestyle="solid"; }

			if(t<0){ color = "red"; }
			else { color = "blue"; }
			linewidth = my::tostring(std::abs(t))+"mm";

			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1),"linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);

			ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+0.2,"\\tiny{"+my::tostring(t)+"}");
		}

		mu = H_(s0,s0);
		if(std::abs(mu)>1e-4){
			if(mu<0){ color = "magenta"; }
			else { color = "cyan"; }
			ps.circle(xy0,std::abs(mu),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
		}

		if(obs_[0].nval()){/*bound energy*/
			corr = obs_[0][i].get_x();
			if(std::abs(corr)>1e-4){
				if(corr<0){ color = "red"; }
				else { color = "blue"; }
				linewidth = my::tostring(std::abs(corr))+"mm";

				ps.line("-",xy0(0)+x_shift,xy0(1),xy1(0)+x_shift,xy1(1),"linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);

				ps.put((xy0(0)+xy1(0))/2.0+x_shift,xy0(1)+0.2,"\\tiny{"+my::tostring(corr).substr(0,5)+"}");
			}
		}

		ps.put(xy0(0),-0.2,"\\tiny{"+my::tostring(s0)+"}");
	}
	if(obs_.size()==2){/*long range correlations*/
		xy0(1) = -1.5;
		xy1(1) = -1.5;
		double rescale(0.75/obs_[1][0].get_x());
		unsigned int n(std::min(4*N_/m_+1,obs_[1].nval()));
		unsigned int idx;
		for(unsigned int i(0);i<n;i++){
			idx = (obs_[1].nval()-n/2+i)%obs_[1].nval();

			corr = obs_[1][idx].get_x()*rescale;
			xy0(0) = i;
			if(std::abs(corr)>1e-4){
				if(i!=idx){
					if(corr<0){ color = "red"; }
					else { color = "blue"; }
				} else { color = "black"; }
				ps.circle(xy0,std::abs(corr),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
			}
		}
	}
	ps.end(true,true,true);
}
/*}*/
