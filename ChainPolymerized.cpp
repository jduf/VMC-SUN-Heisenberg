#include "ChainPolymerized.hpp"

ChainPolymerized::ChainPolymerized(System const& s, Vector<double> const& t):
	System(s),
	Chain<double>(set_spuc(t,N_/m_),"chain-polymerized"),
	t_(t)
{
	if(status_==2){
		init_fermionic();
		filename_ += "-t";
		for(unsigned int i(0);i<t_.size();i++){ filename_ += ((t_(i)>0)?"+":"")+my::tostring(t_(i)); }
		system_info_.text("ChainPolymerized :");
		system_info_.item("Each color has the same Hamiltonian.");
		if(spuc_ != 1){ system_info_.item("Different real hopping terms : "); }
		else { system_info_.item("Uniform real hopping"); }
	}
}

/*{method needed for running*/
void ChainPolymerized::compute_H(){
	H_.set(n_,n_,0);
	unsigned int k(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_:1)*t_(k);
		k = (k+1)%spuc_;
	}
	H_ += H_.transpose();
}

void ChainPolymerized::create(){
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
	if(status_==2 && spuc_!=1){
		/*!Use the eigenvector (k1+k2)/sqrt(2) which correspond to the
		 * impulsion k1+k2=0.*/
		compute_H();
		diagonalize(false);
		double n1(0);
		double n2(0);
		std::complex<double> tmp1;
		std::complex<double> tmp2;
		unsigned int m(M_(0));
		for(unsigned int i(0);i<n_;i++){
			tmp1 = evec_(i,m) + evec_(i,m-1);//k=k1+k2=0
			tmp2 = evec_(i,m) - evec_(i,m-1);//k=k1-k2=2k1
			evec_(i,m-1)= tmp1;
			evec_(i,m)  = tmp2;
			n1 += my::norm_squared(tmp1);
			n2 += my::norm_squared(tmp2);
		}
		for(unsigned int i(0);i<n_;i++){
			evec_(i,m-1)/= sqrt(n1);
			evec_(i,m)  /= sqrt(n2);
		}
		for(unsigned int c(0);c<N_;c++){
			for(unsigned int i(0);i<n_;i++){
				EVec_[c](i,M_(c)-1) = real(evec_(i,M_(c)-1));
			}
		}
	}
}

void ChainPolymerized::save_param(IOFiles& w) const {
	if(w.is_binary()){
		w<<t_;
		w.add_to_header()->title("t=("+my::tostring(t_)+")",'<');
		w.add_to_header()->add(system_info_.get());
	} else { w<<t_<<" "; }
}

unsigned int ChainPolymerized::set_spuc(Vector<double> const& t, unsigned int const& spuc){
	if(t.size() == spuc){ return spuc; }
	else {
		std::cerr<<__PRETTY_FUNCTION__<<" : invalid t size : "<<t.size()<<std::endl;
		return 1;
	}
}
/*}*/

/*{method needed for checking*/
void ChainPolymerized::check(){
	plot_bond_energy();
	plot_long_range_correlations_and_structure_factor();
}

void ChainPolymerized::plot_bond_energy(){
	IOFiles corr_file(analyse_+path_+dir_+filename_+"-corr.dat",true,false);
	corr_file<<"%(2i+1)/2 corr(i,i+1) dx conv(0|1) #conv mean(0|1)"<<IOFiles::endl;

	Vector<double> poly_e(N_/m_,0);
	for(unsigned int i(0);i<obs_[0].nval();i++){
		corr_file<<i+0.5<<" "<<obs_[0][i]<<IOFiles::endl;
		poly_e(i%(N_/m_)) += obs_[0][i].get_x();
	}
	poly_e /= n_*m_/N_;
	poly_e.sort(std::less<double>());

	Gnuplot gp(analyse_+path_+dir_,filename_+"-corr");
	gp+="set key center";
	gp.label("x","site","offset 0,0.5");
	gp.label("y2","$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$");
	gp+="plot '"+filename_+"-corr.dat' u 1:2:3 w errorbars lt 1 lc 7 notitle,\\";
	gp+="     "+my::tostring(poly_e(N_/m_-1)) + " w l lc 3 t 'd-merization="+my::tostring(poly_e(N_/m_-1)-poly_e(N_/m_-2))+"',\\";
	gp+="     "+my::tostring(poly_e(N_/m_-2)) + " w l lc 3 notitle";
	gp.save_file();
	gp.create_image(true,"png");

	if(jd_write_){ jd_write_->write("polymerization strength",poly_e(N_/m_-1)-poly_e(N_/m_-2)); }
}

void ChainPolymerized::display_results(){
	plot_bond_energy();
	plot_long_range_correlations_and_structure_factor();
	if(rst_file_){
		rst_file_->figure(dir_+filename_+"-corr.png",RST::math("E="+my::tostring(obs_[0][0].get_x())+"\\pm"+my::tostring(obs_[0][0].get_dx())),RST::target(dir_+filename_+"-corr.gp")+RST::scale("200"));
		rst_file_->figure(dir_+filename_+"-lr.png","long range correlations",RST::target(dir_+filename_+"-lr.gp")+RST::scale("200"));
		rst_file_->figure(dir_+filename_+"-sf.png","structure factor",RST::target(dir_+filename_+"-sf.gp")+RST::scale("200"));
	}
}
/*}*/

/*{method needed for analysing*/
std::string ChainPolymerized::extract_level_8(){
	rst_file_ = new RSTFile(info_+path_+dir_,filename_);
	std::string basename("../../../../../../../../"+analyse_+path_+dir_+filename_);
	std::string t_string("(");
	for(unsigned int i(0);i<t_.size()-1;i++){
		t_string += my::tostring(t_(i))+",";
	}
	t_string += my::tostring(t_.back())+")";
	std::string title("$N="+my::tostring(N_)+"$ $m="+my::tostring(m_)+"$ $n="+my::tostring(n_)+"$ bc="+my::tostring(bc_)+" $t_{ij}="+t_string+"$");

	(*data_write_)<<t_<<" "<<obs_[0][0]<<IOFiles::endl;
	jd_write_->add_to_header()->title("System's parameters",'-');
	save_param(*jd_write_);
	save(*jd_write_);

	plot_bond_energy();
	rst_file_->figure(basename+"-corr.png","Correlation on links",RST::target(basename+"-corr.gp")+RST::width("1000"));

	plot_long_range_correlations_and_structure_factor();
	rst_file_->figure(basename+"-long-range-corr.png","Long range correlation",RST::target(basename+"-long-range-corr.gp")+RST::width("1000"));
	rst_file_->figure(basename+"-structure-factor.png","Structure factor",RST::target(basename+"-structure-factor.gp")+RST::width("1000"));

	rst_file_->text(read_->get_header());
	rst_file_->save(true,false,true);
	delete rst_file_;
	rst_file_ = NULL;

	return t_string;
}

std::string ChainPolymerized::extract_level_7(){
	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp.title("$N="+my::tostring(N_)+"$ $m="+my::tostring(m_)+"$ $n="+my::tostring(n_)+"$ bc=$"+my::tostring(bc_)+"$");
	if(N_/m_!=4){
		gp.label("x","$t_"+my::tostring(N_/m_)+"$"," offset 0,1");
		gp.label("y2","$\\dfrac{E}{n}$","rotate by 0");
		gp.range("x","0.0","");
		gp+="f(x) = "+std::string(my::are_equal(t_(N_/m_-1),0)?"a+b*x**c":"a+b*(x-c)*(x-c)");
		gp+="a="+my::tostring(obs_[0][0].get_x());
		gp+="b=1";
		gp+="c=1";
		gp+="set fit quiet";
		gp+="fit f(x) '"+filename_+".dat' u "+my::tostring(N_/m_)+":($" +my::tostring(N_/m_+5)+"==0?$"                              +my::tostring(N_/m_+1)+":1/0):"+my::tostring(N_/m_+2)+" zerror via a,b,c";
		gp+="plot '"+filename_+".dat' u "    +my::tostring(N_/m_)+":($" +my::tostring(N_/m_+5)+"==0?$"                              +my::tostring(N_/m_+1)+":1/0):"+my::tostring(N_/m_+2)+" lc 7 w e notitle";
		gp+="     f(x) lc 7 lw 0.5 "+std::string(my::are_equal(t_(N_/m_-1),0)?"notitle":"t sprintf('min %3.4f',c)");
	} else {
		gp.label("x","$t_2$");
		gp.label("y","$t_4$");
		gp.label("z","$\\dfrac{E}{n}$","rotate by 0");
		gp+="set xyplane 0";
		gp+="splot  '"+filename_+".dat' u 2:4:($9==0?$5:1/0) lc 7 notitle";
	}

	gp.save_file();
	gp.create_image(true,"png");

	jd_write_->add_to_header()->title("System's parameters",'-');
	save_param(*jd_write_);
	save(*jd_write_);
	jd_write_->write("polymerization strength",read_->read<double>());
	jd_write_->write("critical exponents",read_->read<Vector<double> >());

	return filename_;
}
/*}*/
