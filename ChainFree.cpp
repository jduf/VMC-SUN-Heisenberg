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
void ChainFree::compute_H(){
	H_.set(n_,n_,0);

	unsigned int s0(0);
	unsigned int s1(0);
	unsigned int k(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(k);
		H_(s0,s0) = mu_(k)/2.0;
		k = (k+1)%spuc_;
	}
	H_ += H_.transpose();
}

void ChainFree::create(){
	compute_H();
	diagonalize(true);
	if(status_==1){
		for(unsigned int c(0);c<N_;c++){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		}
	}
}

void ChainFree::save_param(IOFiles& w) const {
	///*{Description*/
	//std::string t_string("");
	//for(unsigned int i(0);i<t_.size()-1;i++){ t_string += my::tostring(t_(i))+","; }
	//t_string += my::tostring(t_.back());
	//w.write("t ("+t_string+")",t_);
//
	//std::string mu_string("");
	//for(unsigned int i(0);i<mu_.size()-1;i++){ mu_string += my::tostring(mu_(i))+","; }
	//mu_string += my::tostring(mu_.back());
	//w.write("mu ("+mu_string+")",mu_);
	///*}*/
	Vector<double> param(t_.size()+mu_.size());
	for(unsigned int i(0);i<t_.size();i++){ param(i) = t_(i); }
	for(unsigned int i(0);i<mu_.size();i++){ param(i+t_.size()) = mu_(i); }
	
	w.add_header()->title("param (t,mu)",'<');
	w<<param;
	GenericSystem<double>::save_param(w);
}

unsigned int ChainFree::set_spuc(Vector<double> const& t, Vector<double> const& mu, unsigned int const& spuc){
	if(t.size() == mu.size() && mu.size()%spuc==0 && !my::are_equal(t,Vector<double>(spuc,1.0))){ return spuc; }
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
	compute_H();

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

			if(t>0){ color = "blue"; }
			else   { color = "red"; }
			linewidth = my::tostring(std::abs(t))+"mm";

			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1),"linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
		ps.put((xy0(0)+xy1(0))/2.0,xy0(1)+0.2,"\\tiny{"+my::tostring(t)+"}");

		mu = H_(s0,s0);
		if(std::abs(mu)>1e-4){
			if(mu<0){ color = "magenta"; }
			else { color = "cyan"; }
			ps.circle(xy0,std::abs(mu),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
		}
		ps.put(xy0(0),-0.2,"\\tiny{"+my::tostring(mu)+"}");

		if(obs_[0].nval()){/*bound energy*/
			corr = obs_[0][i].get_x();
			if(std::abs(corr)>1e-4){
				if(corr>0){ color = "blue"; }
				else      { color = "red"; }
				linewidth = my::tostring(std::abs(corr))+"mm";

				ps.line("-",xy0(0)+x_shift,xy0(1),xy1(0)+x_shift,xy1(1),"linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);

				ps.put((xy0(0)+xy1(0))/2.0+x_shift,xy0(1)+0.2,"\\tiny{"+my::tostring(corr).substr(0,5)+"}");
			}
		}
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
					if(corr>0){ color = "blue"; }
					else      { color = "red"; }
				} else { color = "black"; }
				ps.circle(xy0,std::abs(corr),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
			}
		}
	}
	ps.end(true,true,true);
}

void ChainFree::display_results(){
	energy_bound();
	long_range_correlation_and_structure_factor();
	if(rst_file_){
		std::string relative_path(analyse_+path_+dir_);
		unsigned int a(std::count(relative_path.begin()+1,relative_path.end(),'/')-1);
		for(unsigned int i(0);i<a;i++){ relative_path = "../"+relative_path; }

		std::string title("J=(");
		for(unsigned int i(0);i<spuc_-1;i++){ title += my::tostring(J_(i)) + ","; }
		title += my::tostring(J_(spuc_-1)) + "), t=(";
		for(unsigned int i(0);i<t_.size()-1;i++){ title += my::tostring(t_(i)) + ","; }
		title += my::tostring(t_.back()) + "), "+RST::math("\\mu")+"=(";
		for(unsigned int i(0);i<mu_.size()-1;i++){ title += my::tostring(mu_(i)) + ","; }
		title += my::tostring(mu_.back()) + ")";
		rst_file_->title(title,'-');

		rst_file_->figure(dir_+filename_+"-pstricks.png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(relative_path+filename_+"-pstricks.pdf")+RST::scale("200"));
		rst_file_->figure(relative_path+filename_+"-lr.png","long range correlations",RST::target(relative_path+filename_+"-lr.gp")+RST::scale("200"));
		rst_file_->figure(relative_path+filename_+"-sf.png","structure factor",RST::target(relative_path+filename_+"-sf.gp")+RST::scale("200"));
	}
}
/*}*/
