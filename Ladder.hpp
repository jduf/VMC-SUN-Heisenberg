#ifndef DEF_LADDER
#define DEF_LADDER

#include "System1D.hpp"

/*{Description*/
/*!
 * This is our ladder.
 *
 * 		1___3___5___7__...__n-1
 * 		|   |   |   |        |
 * 		0___2___4___6__...__n-2
 */
/*}*/
template<typename Type>
class Ladder: public System1D<Type>{
	public:
		/*{Description*/
		/*!Constructor of a system with spuc sites per unit cell. Calls the
		 * GenericSystem<Type>(3,filename), to construct a system with 3 links
		 * per sites */
		/*}*/
		Ladder(unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Ladder()=0;

	protected:
		/*!Implement the pure virtual method in GenericSystem*/
		void long_range_correlations_obs();
		/*!Returns the neighbours of site i*/
		Matrix<int> get_neighbourg(unsigned int const& i) const;
		/*Draw the lattice inside a PSTricks file*/
		void draw_lattice(bool const& only_unit_cell, bool const& silent, bool const& create_image, std::string param, std::string title);
		/*!Computes and writes the flux per plaquette in the PSTricks file*/
		std::string flux_per_plaquette(unsigned int const& s0, unsigned int const& s1) const;
		/*!To get the data for a plot of E(theta) for all N,m,n,bc*/
		std::string extract_level_6();
		/*!Given N and m, save the best simulation in a text file for any n*/
		std::string extract_level_3();
};

/*{constructor*/
template<typename Type>
Ladder<Type>::Ladder(unsigned int const& spuc, std::string const& filename):
	System1D<Type>(spuc,3,filename)
{
	/*!(*this) has been created via System(System const& s), this->J_ and
	 * this->links_ are already a copy of s. therefore if s has correctly
	 * defined J_ and links_, there is no need to recompute them for (*this).
	 * if s has undefined links_ and if J_ is of size 2, then this->links_ and
	 * this_->J_ should be computed*/
	if(this->status_==3){
		/*!create the links if necessary*/
		if(this->ref_(4)==2){
			Vector<unsigned int> tmp(2);
			tmp(0) = 2;
			tmp(1) = 1;
			this->energy_obs(tmp);
		} else {
			this->ref_(4) = 0;
			this->status_ = 2;
		}
		

		/*!sets the bond energy if it has not been set yet*/
		if(this->obs_[0].nlinks() != this->J_.size() && this->J_.size() == 2){
			Vector<double> tmp(this->J_);
			this->J_.set(this->obs_[0].nlinks());
			for (unsigned int i(0); i<this->J_.size();i++){
				if(this->obs_[0](i,3)){ this->J_(i) = tmp(1); } //rungs (J⊥) -> sin(theta)
				else { this->J_(i) = tmp(0); } //legs  (J‖) -> cos(theta)
			}
		}

		/*!fix the names for the bond energy*/
		if(this->J_.size()==this->obs_[0].nlinks()){
			std::string tmp("theta"+my::tostring(acos(this->J_(0))));
			this->filename_.replace(this->filename_.find("Juniform"),8,tmp);
			this->path_.replace(this->path_.find("Juniform"),8,tmp);
		} else { std::cerr<<__PRETTY_FUNCTION__<<" : J_ has an incoherent size"<<std::endl; }
	} else {
		/*!if the ladder has a spuc_ equal to one, the creation is impossible
		 * but the construction should be silent to have a nice display when
		 * used with Analyse* */
		if(this->spuc_!=1){ std::cerr<<__PRETTY_FUNCTION__<<" creation is problematic"<<std::endl; }
	}
}

template<typename Type>
Ladder<Type>::~Ladder() = default;
/*}*/

/*{protected methods*/
template<typename Type>
void Ladder<Type>::long_range_correlations_obs(){
	unsigned int m(this->n_/2);
	unsigned int nval(this->n_/2);
	unsigned int nlinks(m*nval);
	unsigned int idx(this->obs_.size());
	this->obs_.push_back(Observable("S_10*S1i",2,nval,nlinks));
	this->obs_.push_back(Observable("S_10*S2i",2,nval,nlinks));
	this->obs_.push_back(Observable("S_20*S1i",2,nval,nlinks));
	this->obs_.push_back(Observable("S_20*S2i",2,nval,nlinks));
	for(unsigned int i(0);i<m;i++){
		for(unsigned int j(0);j<nval;j++){
			/*obs_[1]=S_10*S_1i*/
			this->obs_[idx](i*nval+j,0) = 2*i;
			this->obs_[idx](i*nval+j,1) = (2*(i+j))%this->n_;
			this->obs_[idx](i*nval+j,2) = j;
			/*obs_[2]=S_10*S_2i*/
			this->obs_[idx+1](i*nval+j,0) = 2*i;
			this->obs_[idx+1](i*nval+j,1) = (2*(i+j)+1)%this->n_;
			this->obs_[idx+1](i*nval+j,2) = j;
			/*obs_[3]=S_20*S_1i*/
			this->obs_[idx+2](i*nval+j,0) = 2*i+1;
			this->obs_[idx+2](i*nval+j,1) = (2*(i+j))%this->n_;
			this->obs_[idx+2](i*nval+j,2) = j;
			/*obs_[4]=S_20*S_2i*/
			this->obs_[idx+3](i*nval+j,0) = 2*i+1;
			this->obs_[idx+3](i*nval+j,1) = (2*(i+j)+1)%this->n_;
			this->obs_[idx+3](i*nval+j,2) = j;
		}
	}
}

template<typename Type>
Matrix<int> Ladder<Type>::get_neighbourg(unsigned int const& i) const {
	Matrix<int> nb(this->z_,3);
	if(i%2){//!odd number => upper part. 0:right, 1:down, 2:left
		if(i+1 != this->n_){
			nb(0,0) = i+2;
			nb(0,1) = 0;
			nb(0,2) = 0;
		} else {
			nb(0,0) = 1;
			nb(0,1) = 0;
			nb(0,2) = 1;
		}
		nb(1,0) = i-1;
		nb(1,1) = 1;
		nb(1,2) = 0;
		if(i != 1){
			nb(2,0) = i-2;
			nb(2,1) = 2;
			nb(2,2) = 0;
		} else {
			nb(2,0) = this->n_-1;
			nb(2,1) = 2;
			nb(2,2) = 1;
		}
	} else {//!even number => lower part. 0:right, 1:up, 2:left
		if(i+2 != this->n_){
			nb(0,0) = i+2;
			nb(0,1) = 0;
			nb(0,2) = 0;
		} else{
			nb(0,0) = 0;
			nb(0,1) = 0;
			nb(0,2) = 1;
		}
		nb(1,0) = i+1;
		nb(1,1) = 1;
		nb(1,2) = 0;
		if(i){
			nb(2,0) = i-2;
			nb(2,1) = 2;
			nb(2,2) = 0;
		} else {
			nb(2,0) = this->n_-2;
			nb(2,1) = 2;
			nb(2,2) = 1;
		}
	}
	return nb;
}

template<typename Type>
void Ladder<Type>::draw_lattice(bool const& only_unit_cell, bool const& silent, bool const& create_image, std::string param, std::string title){
	Matrix<int> links(this->obs_[0].get_links());
	Vector<unsigned int> o(6,0);
	double max_bond_energy(0);
	unsigned int o_index_tmp(2);
	for(unsigned int i(1);i<this->obs_.size();i++){
		switch(this->obs_[i].get_type()){
			case 1:
				{
					o(0)=i;
					for(unsigned int j(0);j<this->obs_[i].nval();j++){
						if(max_bond_energy < std::abs(this->obs_[i][j].get_x()/(this->m_*this->m_))){
							max_bond_energy = std::abs(this->obs_[i][j].get_x()/(this->m_*this->m_));
						}
					}
				}break;//bond energy
			case 2:{ o(o_index_tmp++)=i; }break;//long range correlation
			case 3:{ o(1)=i; }break;//color occupation
		}
	}

	std::string arrow("-");
	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	unsigned int s0;
	unsigned int s1;
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	Type t;
	double mu;
	double bond_energy;
	Vector<double> shift(2,0.0);
	Vector<unsigned int> loop(4);
	loop(0) = 0;
	loop(1) = 1;
	loop(2) = 2;
	loop(3) = 3;
	Matrix<double> uc(4,2);
	uc(0,0) =-0.5;
	uc(0,1) =-0.5;
	uc(1,0) = this->spuc_/2.0-0.5;
	uc(1,1) =-0.5;
	uc(2,0) = this->spuc_/2.0-0.5;
	uc(2,1) = 1.5;
	uc(3,0) =-0.5;
	uc(3,1) = 1.5;
	PSTricks ps(this->get_info_path(),this->filename_);
	ps.add("\\newcommand{\\wbg}[1]{\\setlength{\\fboxsep}{ 1pt}\\colorbox{white}{\\tiny{#1}}}");
	ps.begin(-1,-5,this->n_/1.5,2,this->filename_);
	ps.polygon(uc,"linecolor=black,linestyle=dashed");
	unsigned int link_type;
	if(only_unit_cell){
		shift(0) = uc(1,0)-uc(0,0)+0.5;
		for(unsigned int i(0);i<links.row();i++){
			s0 = links(i,0);
			s1 = links(i,1);
			if(!(s0%2)){
				if(!(s1%2)){ link_type = 0; } // lower horizontal link
				else{ link_type = 1; } // vertical link
			} else { link_type = 2; } // upper horizontal link
			/*!shift xy 0 to avoid latex bug for value too close to 0*/
			xy0(0) = s0/2+0.0001;
			xy0(1) = s0%2+0.0001;
			xy1(0) = s1/2+0.0001;
			xy1(1) = s1%2+0.0001;

			if(links(i,4)){
				linestyle = "dashed";
				xy1(0) = xy0(0)+1.0;
			} else { linestyle = "solid"; }

			if(my::in_polygon(uc.row(),uc.ptr(),uc.ptr()+uc.row(),xy0(0),xy0(1))){
				/*Hopping amplitude, chemical potential and fluxes*/
				t = this->H_(s0,s1);
				if(std::abs(t)>1e-4){
					linewidth = my::tostring(std::abs(t))+"mm";
					if(my::real(t)>0){ color = "blue"; }
					else             { color = "red"; }
					if(my::are_equal(my::imag(t),0)){ arrow = "-"; }
					else {
						if(my::imag(-t)>0){ arrow = "->"; }
						else              { arrow = "<-"; }
					}
					ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle=solid");
					ps.put((xy0(0)+xy1(0))/2.0,(xy0(1)+xy1(1))/2.0, "\\wbg{"+my::tostring(my::round_nearest(std::abs(t),1000))+"}");
				}
				if(link_type!=1){
					mu = my::real(this->H_(s0,s0))/M_PI;
					if(std::abs(mu)>1e-4){
						if(mu>0){ color = "cyan"; }
						else    { color = "magenta"; }
						ps.circle(xy0,sqrt(std::abs(mu)),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
					}
					if(link_type==0){ ps.put(xy0(0)+0.5,xy0(1)+0.5,flux_per_plaquette(s0,s1)); }
				}

				/*Bond energy and color occupation*/
				xy0 += shift;
				xy1 += shift;
				if(o(0)){
					bond_energy = this->obs_[o(0)][links(i,2)].get_x()/(this->m_*this->m_);
					linewidth = my::tostring(std::abs(bond_energy))+"mm";
					if(std::abs(bond_energy)>1e-4){
						if(bond_energy>0){ color = "blue"; }
						else             { color = "red"; }
						ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle=solid");
					}
					ps.put((xy0(0)+xy1(0))/2.0,(xy0(1)+xy1(1))/2.0, "\\wbg{"+my::tostring(my::round_nearest(std::abs(bond_energy)/max_bond_energy,100))+"}");
				}
				if(link_type!=1 && o(1)){
					xy0 += shift;
					Vector<double> p(this->N_);
					for(unsigned int j(0);j<this->N_;j++){ p(j) = this->obs_[o(1)][j+this->N_*links(i,5)].get_x(); }
					ps.pie(xy0(0),xy0(1),p,0.2,"chartColor=color");
				}
			}
		}
	} else {
		for(unsigned int i(0);i<links.row();i++){
			s0 = links(i,0);
			s1 = links(i,1);
			if(!(s0%2)){
				if(!(s1%2)){ link_type = 0; } // lower horizontal link
				else{ link_type = 1; } // vertical link
			} else { link_type = 2; } // upper horizontal link
			/*!shift xy 0 to avoid latex bug for value too close to 0*/
			xy0(0) = s0/2+0.0001;
			xy0(1) = s0%2+0.0001-2;
			xy1(0) = s1/2+0.0001;
			xy1(1) = s1%2+0.0001-2;

			if(links(i,4)){
				linestyle = "dashed";
				xy1(0) = xy0(0)+1.2;
				xy1(0) = xy0(0)+1.0;
			} else { linestyle = "solid"; }

			/*Draws only the lattice, shows links and bc and indices*/
			linewidth = my::tostring(this->J_(i)/this->J_.max()*1)+"pt";
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1),"linewidth="+linewidth+",linecolor=black,linestyle="+linestyle);
			if(link_type==1){
				ps.put(xy0(0)-0.15,xy0(1)-0.15,"\\tiny{"+my::tostring(s0)+"}");
				ps.put(xy0(0)+0.15,xy0(1)-0.15,"\\textcolor{green}{\\tiny{"+my::tostring(links(i,5))+"}}");
				ps.put(xy0(0)-0.15,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
				ps.put(xy0(0)+0.15,xy1(1)+0.15,"\\textcolor{green}{\\tiny{"+my::tostring(links(i,6))+"}}");
			}

			/*Bond energy and color occupation*/
			xy0(1) -= 2.0;
			xy1(1) -= 2.0;
			if(o(0)){
				bond_energy = this->obs_[o(0)][links(i,2)].get_x()/(this->m_*this->m_);
				linewidth = my::tostring(std::abs(bond_energy))+"mm";
				if(std::abs(bond_energy)>1e-4){
					if(bond_energy>0){ color = "blue"; }
					else             { color = "red"; }
					ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle=solid");
				}
			}
			if(link_type!=1 && o(1)){
				Vector<double> p(this->N_);
				for(unsigned int j(0);j<this->N_;j++){ p(j) = this->obs_[o(1)][j+this->N_*links(i,5)].get_x(); }
				ps.pie(xy0(0),xy0(1),p,0.2,"chartColor=color");
			}

			/*Hopping amplitude, chemical potential and fluxes*/
			xy0(1) += 4.0;
			xy1(1) += 4.0;
			t = this->H_(s0,s1);
			if(std::abs(t)>1e-4){
				linewidth = my::tostring(std::abs(t))+"mm";
				if(my::real(t)>0){ color = "blue"; }
				else             { color = "red"; }
				if(my::are_equal(my::imag(t),0)){ arrow = "-"; }
				else {
					if(my::imag(-t)>0){ arrow = "->"; }
					else              { arrow = "<-"; }
				}
				ps.line(arrow,xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
			}
			if(link_type!=1){
				mu = my::real(this->H_(s0,s0))/M_PI;
				if(std::abs(mu)>1e-4){
					if(mu>0){ color = "cyan"; }
					else    { color = "magenta"; }
					ps.circle(xy0,sqrt(std::abs(mu)),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
				}
				if(link_type==0){ ps.put(xy0(0)+0.5,xy0(1)+0.5,flux_per_plaquette(s0,s1)); }
			}
		}
	}
	if(o(2) && o(3)){/*long range correlations*/
		double corr;
		double rescale(0.75/this->obs_[o(2)][0].get_x());
		unsigned int n(std::min(4*this->N_/this->m_+1,this->obs_[o(2)].nval()));
		unsigned int idx;
		for(unsigned int i(0);i<n;i++){
			idx = (this->obs_[o(2)].nval()-n/2+i)%this->obs_[o(2)].nval();

			corr = this->obs_[o(2)][idx].get_x()*rescale;
			xy0(0) = i;
			xy0(1) =-2;
			if(std::abs(corr)>1e-4){
				if(i!=idx){
					if(corr>0){ color = "blue"; }
					else      { color = "red"; }
				} else { color = "black"; }
				ps.circle(xy0,std::abs(corr),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
			}

			corr = this->obs_[o(3)][idx].get_x()*rescale;
			xy0(1) =-1;
			if(std::abs(corr)>1e-4){
				if(corr>0){ color = "blue"; }
				else      { color = "red"; }
				ps.circle(xy0,std::abs(corr),"fillstyle=solid,fillcolor="+color+",linecolor="+color);
			}
		}
	}
	ps.end(silent,true,true);

	/*!long range correlations*/
	if(o(2) && o(3) && o(4) && o(5)){
		unsigned int llr(this->obs_[o(2)].nval());
		Vector<std::complex<double> > Ck_ll(llr,0.0);
		Vector<std::complex<double> > Ck_llp(llr,0.0);
		Vector<std::complex<double> > Ck_m(llr,0.0);
		Vector<std::complex<double> > Ck_p(llr,0.0);
		double dk(2.0*M_PI/llr);

		for(unsigned int i(0);i<llr;i++){
			for(unsigned int k(0);k<llr;k++){
				Ck_ll(k) += std::polar(this->obs_[o(2)][i].get_x()+this->obs_[o(5)][i].get_x(),dk*k*i);
				Ck_llp(k)+= std::polar(this->obs_[o(3)][i].get_x()+this->obs_[o(4)][i].get_x(),dk*k*i);
				Ck_m(k)  += std::polar(this->obs_[o(2)][i].get_x()-this->obs_[o(3)][i].get_x()-this->obs_[o(4)][i].get_x()+this->obs_[o(5)][i].get_x(),dk*k*i);
				Ck_p(k)  += std::polar(this->obs_[o(2)][i].get_x()+this->obs_[o(3)][i].get_x()+this->obs_[o(4)][i].get_x()+this->obs_[o(5)][i].get_x(),dk*k*i);
			}
		}
		double sum_ll(Ck_ll.sum().real());
		double sum_llp(Ck_llp.sum().real());
		double sum_m(Ck_m.sum().real());
		double sum_p(Ck_p.sum().real());
		Ck_ll /= 0.5*dk*sum_ll;
		Ck_llp/= 0.5*dk*sum_llp;
		Ck_m  /= 0.5*dk*sum_m;
		Ck_p  /= 0.5*dk*sum_p;

		IOFiles file_co(this->analyse_+this->path_+this->dir_+this->filename_+"-lrc.dat",true,false);
		IOFiles file_sf(this->analyse_+this->path_+this->dir_+this->filename_+"-sf.dat",true,false);
		for(unsigned int l(0);l<llr;l++){
			file_co<<l<<" "<<this->obs_[o(2)][l]<<" "<<this->obs_[o(3)][l]<<" "<<this->obs_[o(4)][l]<<" "<<this->obs_[o(5)][l]<<IOFiles::endl;
			file_sf<<dk*l<<" "<<Ck_m(l).real()<<" "<<Ck_m(l).imag()<<" "<<sum_m<<" "<<Ck_p(l).real()<<" "<<Ck_p(l).imag()<<" "<<sum_p<<" "<<Ck_ll(l).real()<<" "<<Ck_ll(l).imag()<<" "<<sum_ll<<" "<<Ck_llp(l).real()<<" "<<Ck_llp(l).imag()<<" "<<sum_llp<<IOFiles::endl;
		}

		/*!inter (intra) correlation and structure factor*/
		/*{*/
		Gnuplot gp_in(this->analyse_+this->path_+this->dir_,this->filename_+"-lr");
		gp_in+="data_co='"+this->filename_+"-lrc.dat'";
		gp_in+="data_sf='"+this->filename_+"-sf.dat'";
		gp_in.multiplot();
		/*{correlations*/
		gp_in.range("x","0",llr/2);

		gp_in.tics("x");
		gp_in.margin("0.1","0.5","0.9","0.5");
		gp_in+="plot data_co u 1:($2+$14):($3+$18) w errorbars lt 1 lc 6 t '$\\tilde{C}^{ll}(k_x)$'";

		gp_in.margin("0.1","0.5","0.5","0.1");
		gp_in.tics("x","");
		gp_in+="plot data_co u 1:($6+$10):($7+$11)  w errorbars lt 1 lc 7 t '$\\tilde{C}^{llp}(k_x)$'";
		/*}*/
		/*{structure factor*/
		gp_in.range("x","0","pi");

		gp_in.key("left");
		gp_in.tics("x");
		gp_in.tics("y");
		gp_in.tics("x2","('' pi/3, '' pi/2, '' 2*pi/3, '' pi) mirror");
		gp_in.tics("y2","mirror");
		gp_in.margin("0.5","0.9","0.9","0.5");
		gp_in+="plot data_sf u 1:8 axes x1y2 lt "+std::string(sum_ll<0?"1":"2")+" lc 6 notitle,\\";
		gp_in+="     data_sf u 1:9 axes x1y2 lt "+std::string(sum_ll<0?"1":"2")+" lc 6 notitle";

		gp_in.margin("0.5","0.9","0.5","0.1");
		gp_in.tics("x","('$\\pi/3$' pi/3, '$\\pi/2$' pi/2, '$2\\pi/3$' 2*pi/3, '$\\pi$' pi)");
		gp_in+="plot data_sf u 1:11 axes x1y2 lt "+std::string(sum_llp<0?"1":"2")+" lc 7 notitle,\\";
		gp_in+="     data_sf u 1:12 axes x1y2 lt "+std::string(sum_llp<0?"1":"2")+" lc 7 notitle";
		/*}*/
		gp_in.save_file();
		if(create_image){ gp_in.create_image(silent,"png"); }
		/*}*/
		/*!(anti)symmetric correlations and structure factors*/
		/*{*/
		Gnuplot gp_sas(this->analyse_+this->path_+this->dir_,this->filename_+"-as");
		gp_sas+="data_co='"+this->filename_+"-lrc.dat'";
		gp_sas+="data_sf='"+this->filename_+"-sf.dat'";
		gp_sas.multiplot();
		/*{correlations*/
		gp_sas.range("x","0",llr/2);

		gp_sas.tics("x");
		gp_sas.margin("0.1","0.5","0.9","0.5");
		gp_sas+="plot data_co u 1:($2-$6-$10+$14):($3+$7+$11+$15) w errorbars lt 1 lc 6 t '$-$'";

		gp_sas.margin("0.1","0.5","0.5","0.1");
		gp_sas.tics("x","");
		gp_sas+="plot data_co u 1:($2+$6+$10+$14):($3+$7+$11+$15)  w errorbars lt 1 lc 7 t '$+$'";
		/*}*/
		/*{structure factor*/
		gp_sas.range("x","0","pi");

		gp_sas.key("left");
		gp_sas.tics("x");
		gp_sas.tics("y");
		gp_sas.tics("x2","('' pi/3, '' pi/2, '' 2*pi/3, '' pi)");
		gp_sas.tics("y2","");
		gp_sas.margin("0.5","0.9","0.9","0.5");
		gp_sas+="plot data_sf u 1:2 axes x1y2 lt "+std::string(sum_m<0?"1":"2")+" lc 6 notitle,\\";
		gp_sas+="     data_sf u 1:3 axes x1y2 lt "+std::string(sum_m<0?"1":"2")+" lc 6 notitle";

		gp_sas.margin("0.5","0.9","0.5","0.1");
		gp_sas.tics("x","('$\\pi/3$' pi/3, '$\\pi/2$' pi/2, '$2\\pi/3$' 2*pi/3, '$\\pi$' pi)");
		gp_sas+="plot data_sf u 1:5 axes x1y2 lt "+std::string(sum_p<0?"1":"2")+" lc 7 notitle,\\";
		gp_sas+="     data_sf u 1:6 axes x1y2 lt "+std::string(sum_p<0?"1":"2")+" lc 7 notitle";
		/*}*/
		gp_sas.save_file();
		if(create_image){ gp_sas.create_image(silent,"png"); }
		/*}*/
	}

	if(this->rst_file_){
		title = RST::math("\\theta=")+my::tostring(acos(this->J_(0))) + " : " + title;
		param += " -d:theta " + my::tostring(acos(this->J_(0)));
		this->rst_file_set_default_info(param,title,"theta"+my::tostring(acos(this->J_(0))));

		if(o(2) && o(3) && o(4) && o(5)){
			std::string path(this->analyse_+this->path_+this->dir_);
			std::string rp("");
			int p(0);
			if( (p=path.find('/',p)) == 0){ rp = path; }
			else {
				while( (p=path.find('/',p)+1) ){
					if(path[p-1] != '.' && path[p-1] != '/'){ rp += "../"; }
					else if( p>1 && path[p-1] == '.' && path[p-2] != '/'){ rp += "../"; }
				}
			}
			this->rst_file_->figure(rp+this->filename_+"-lr.png","long range correlations",     RST::target(rp+this->filename_+"-lr.gp")+RST::width("800"));
			this->rst_file_->figure(rp+this->filename_+"-as.png","(anti)symmetric correlations",RST::target(rp+this->filename_+"-as.gp")+RST::width("800"));
		}
	}
}

template<typename Type>
std::string Ladder<Type>::flux_per_plaquette(unsigned int const& s0, unsigned int const& s1) const {
	double flux(0.0);
	flux+= std::arg(-this->H_(s0,s1));
	flux+= std::arg(-this->H_(s1,s1+1));
	flux+= std::arg(-this->H_(s1+1,s0+1));
	flux+= std::arg(-this->H_(s0+1,s0));
	flux /= M_PI;
	if(flux>2.0){ flux -= 2.0; }
	if(flux<-2.0){ flux += 2.0; }

	double sign;
	unsigned long long a;
	unsigned long long b;
	if(!my::are_equal(flux,0.0,this->eq_prec_,this->eq_prec_)){
		if(my::to_fraction(flux,a,b,sign) && b!=1){ return std::string(sign<0?"-":"")+"$\\frac{"+(a==1?"":my::tostring(a))+"\\pi}{"+my::tostring(b)+"}$"; }
		else if((unsigned int)(my::chop(flux))%2){ return "$\\pi$"; }
	}
	return "";
}

template<typename Type>
std::string Ladder<Type>::extract_level_6(){
	(*this->data_write_)<<this->N_<<" "<<this->m_<<" "<<this->bc_<<" "<<this->n_<<" "<<acos(this->J_(0))<<" "<<this->obs_[0][0]<<" "<<this->ref_<<IOFiles::endl;

	return this->filename_;
}

template<typename Type>
std::string Ladder<Type>::extract_level_3(){
	(*this->read_)>>this->obs_[0][0];
	(*this->data_write_)<<this->N_<<" "<<this->m_<<" "<<this->bc_<<" "<<this->n_<<" "<<this->obs_[0][0]<<IOFiles::endl;

	return this->filename_;
}
/*}*/
#endif
