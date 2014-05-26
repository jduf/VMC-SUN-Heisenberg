#include "ChainPolymerized.hpp"

ChainPolymerized::ChainPolymerized(unsigned int const& N, unsigned int const& n, unsigned int const& m, int const& bc, Vector<unsigned int> const& ref):
	Chain<double>(N,n,m,bc,ref,"chain-polymerized")
{
	rst_.text("Spin chain, with different real hopping term.");
	rst_.text("For N colors and m particules per sites, every");
	rst_.text("N/m, there is a weaker bound, namely t-delta");
	rst_.text("instead of t+delta. (t=1,delta>0)");
}

ChainPolymerized::~ChainPolymerized(){}

void ChainPolymerized::create(double const& delta, unsigned int const& type){
	std::cout<<"chainpolymerized create type="<<type<<std::endl;
	delta_=delta;
	T_.set(n_,n_,0);
	EVec_.set(n_*N_,M_,0);

	compute_T();
	diagonalize_T('S');
	for(unsigned int spin(0);spin<N_;spin++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_;j++){
				EVec_(i+spin*n_,j) = T_(i,j);
			}
		}
	}
}

void ChainPolymerized::compute_T(){
	/*!If t<0, delta<0 otherwise no polymerization occurs
	 * If t>0, delta>0 otherwise no polymerization occurs */
	double t(1.0);
	Matrix<int> nb;
	for(unsigned int i(0); i < n_; i += a_){
		for(unsigned int j(0); j<a_-1; j++){
			nb = get_neighbourg(i+j);
			T_(i+j,nb(0,0)) = t+delta_;
		}
		nb = get_neighbourg(i+a_-1);
		T_(i+a_-1,nb(0,0)) = nb(0,1)*(t-delta_);
	}
	T_ += T_.transpose();
}

void ChainPolymerized::compute_P(Matrix<double>& P){
	P.set(n_,n_);
	P(n_ -1,0) = bc_;
	for(unsigned int i(0); i< n_-1; i++){
		P(i,i+1) = 1.0;
	}
}

void ChainPolymerized::save(IOFiles& w) const{
	GenericSystem<double>::save(w);
	w("delta (t+-delta)",delta_);
}

void ChainPolymerized::check(){
	delta_=0.1;
	compute_T();
	std::cout<<T_<<std::endl;
}

void ChainPolymerized::plot_corr(std::string const& path, std::string const& filename, unsigned int const& nruns){
	Gnuplot gp(path,filename+"-corr");
	gp+="set xlabel 'site' offset 0,0.5";
	gp+="set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$' offset 1";
	gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc="+tostring(bc_)+" $\\delta="+tostring(delta_)+"$'";
	gp+="plot for [IDX=0:"+tostring(nruns-1)+"] '"+filename+"-corr.dat' i IDX u 1:($4==1?$2:1/0):3 w errorbars lt 1 lc 3 ps 0 notitle,\\";
	gp+="     for [IDX=0:"+tostring(nruns-1)+"] '"+filename+"-corr.dat' i IDX u 1:($4==0?$2:1/0):3 w errorbars lt 1 lc 1 ps 0 notitle,\\";
	gp+="                   '"+filename+"-corr.dat' i " + tostring(nruns) + " u 1:2:3 w errorbars lt 1 lc 2 lw 2 notitle";
	gp.save_file();
}

void ChainPolymerized::plot_long_range_corr(std::string const& path, std::string const& filename, unsigned int const& nruns){
	unsigned int length(long_range_corr_.size());
	if(length>0){
		Gnuplot gp(path,filename+"-long-range-corr");
		gp+="stats '"+filename+"-long-range-corr.dat' nooutput";
		gp.xrange(0,length+1);
		gp.yrange("1.1*STATS_min_y","1.1*STATS_max_y");
		gp+="set xlabel '$\\|i-j\\|$' offset 0,0.5";
		gp+="set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(j)>$' offset 1";
		gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc="+tostring(bc_)+" $\\delta="+tostring(delta_)+"$'";
		gp+="set key right bottom";
		gp+="a=1.0";
		gp+="b=1.0";
		gp+="eta=1.0";
		gp+="m="+tostring(m_)+".0";
		gp+="N="+tostring(N_)+".0";
		gp+="f(x) = a/(x*x) + b*cos(2.0*pi*x*m/N)/(x**eta)";
		gp+="set fit quiet";
		switch(N_/m_){
			case 2:{ gp+="fit [3:"+tostring(length)+"] f(x) '"+filename+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; } break;
			case 3:{
					   switch((length + 1) % 3){
						   case 0:{ gp+="fit [2:"+tostring(length)+"] f(x) '"+filename+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
						   case 1:{ gp+="fit [5:"+tostring(length)+"] f(x) '"+filename+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
						   case 2:{ gp+="fit [3:"+tostring(length)+"] f(x) '"+filename+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
					   }break;
				   }break;
			default :{ gp+="fit ["+tostring(N_-1)+":] f(x) '"+filename+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
		}
		gp+="plot for [IDX=0:"+tostring(nruns-1)+"] '"+filename+"-long-range-corr.dat' i IDX u 1:($4==1?$2:1/0):3 w errorbars lt 1 lc 3 ps 0 notitle,\\";
		gp+="     for [IDX=0:"+tostring(nruns-1)+"] '"+filename+"-long-range-corr.dat' i IDX u 1:($4==0?$2:1/0):3 w errorbars lt 1 lc 1 ps 0 notitle,\\";
		gp+="                   '"+filename+"-long-range-corr.dat' i "+tostring(nruns)+" u 1:2:3 w errorbars lt 1 lc 2 lw 2 notitle,\\";
		gp+="                   f(x) notitle";
		gp.save_file();
	}
}

void ChainPolymerized::treat_one_sim(IOFiles& read, IOFiles& write, RSTFile& rst, std::string const& path, std::string const& filename){
	unsigned int nruns;
	{
		IOFiles corr_file(path+filename+"-corr.dat",true);
		IOFiles long_range_corr_file(path+filename+"-long-range-corr.dat",true);//should not be delcared when type!=2
		rst.link_figure(path+filename+"-corr.png","Correlation on links",path+filename+"-corr.gp",1000);
		rst.link_figure(path+filename+"-long-range-corr.png","Long range correlation",path+filename+"-long-range-corr.gp",1000);
		read>>nruns;
		for(unsigned int i(0);i<nruns+1;i++){ 
			read>>E_>>corr_>>long_range_corr_;
			for(unsigned int j(0);j<corr_.size();j++){
				corr_file<<j+0.5<<" "<<corr_[j]<<IOFiles::endl;
			}
			for(unsigned int j(0);j<long_range_corr_.size();j++){
				long_range_corr_file<<j+1<<" "<<long_range_corr_[j]<<IOFiles::endl;
			}
			long_range_corr_file<<IOFiles::endl<<IOFiles::endl;
		}
	}
	Vector<double> poly_e(N_/m_,0);
	unsigned int i(0);
	while(i<corr_.size()){
		for(unsigned int j(0);j<N_/m_;j++){
			poly_e(j) += corr_[i].get_x();
			i++;
		}
	}
	poly_e /= n_*m_/N_;
	poly_e.sort();

	//write<<E_<<poly_e(N_/m_-2)-poly_e(N_/m_-1)<<IOFiles::endl;
	write("E",E_);
	write("polymerization strength",poly_e(N_/m_-2)-poly_e(N_/m_-1));
	plot_corr(path,filename,nruns);
	plot_long_range_corr(path,filename,nruns);
}
