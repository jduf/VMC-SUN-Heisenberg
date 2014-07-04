#include "AnalyseParameter.hpp"


void AnalyseParameter::open_files(std::string const& jdfile, std::string const& datafile, Directory const& d){
	if(level_>1){ 
		write_ = new IOFiles(jdfile,true);
		(*write_)("number of jdfiles",d.size());
		write_->add_to_header("\n");
	}
	if(level_==3 || level_==5){
		data_write_ = new IOFiles(datafile,true);
	}
}

void AnalyseParameter::close_files(){
	if(write_){ 
		switch(level_){
			case 5:{ rst_.last().link_figure(analysis_+path_+dir_.substr(0,dir_.size()-1)+".png","E.png",analysis_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000); } break;
			case 3:{ rst_.last().link_figure(analysis_+path_+dir_.substr(0,dir_.size()-1)+".png","d-merization.png",analysis_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000); } break;
		}
		rst_.last().text(write_->get_header());
		delete write_;
		write_ = NULL;
	}
	if(data_write_){
		delete data_write_;
		data_write_ = NULL;
	}
}

void AnalyseParameter::extract_level_5(){
	/*E(param)|n=fixé*/
	IOFiles read(sim_+path_+dir_+filename_+".jdbin",false);
	RSTFile rst(info_+path_+dir_,filename_);

	unsigned int nruns;
	unsigned int tmax;
	double param;

	read>>nruns>>tmax>>ref_>>N_>>m_>>n_>>M_>>bc_>>param;
	IOFiles corr_file(analysis_+path_+dir_+filename_+"-corr.dat",true);
	IOFiles long_range_corr_file(analysis_+path_+dir_+filename_+"-long-range-corr.dat",true);//should not be delcared when type!=2
	data_write_->precision(10);
	(*data_write_)<<"% delta E dE 0|1"<<IOFiles::endl;
	/* the +1 is the averages over all runs */
	Vector<double> poly_e(N_/m_,0);
	for(unsigned int i(0);i<nruns+1;i++){ 
		read>>E_>>corr_>>long_range_corr_;

		if(i<nruns){
			unsigned int k(0);
			while(k<corr_.size()){
				for(unsigned int j(0);j<N_/m_;j++){
					poly_e(j) += corr_[k].get_x();
					k++;
				}
			}
		}

		(*data_write_)<<param<<" "<<E_.get_x()<<" "<<E_.get_dx()<<" "<<(i<nruns?true:false)<<IOFiles::endl;
		for(unsigned int j(0);j<corr_.size();j++){
			corr_file<<j+0.5<<" "<<corr_[j]<<" "<<(i<nruns?true:false)<<IOFiles::endl;
		}
		for(unsigned int j(0);j<long_range_corr_.size();j++){
			long_range_corr_file<<j+1<<" "<<long_range_corr_[j]<<" "<<(i<nruns?true:false)<<IOFiles::endl;
		}

	}
	poly_e /= nruns*n_*m_/N_;
	poly_e.sort(std::less<double>());

	(*write_)("ref (type of wavefunction)",ref_);
	(*write_)("N (N of SU(N))",N_);
	(*write_)("m (# of particles per site)",m_);
	(*write_)("n (# of site)",n_);
	(*write_)("bc (boundary condition)",bc_);
	(*write_)("param",param);
	(*write_)("E",E_);
	(*write_)("polymerization strength",poly_e(N_/m_-1)-poly_e(N_/m_-2));

	/*{*/
	Gnuplot gp(analysis_+path_+dir_,filename_+"-corr");
	gp+="set xlabel 'site' offset 0,0.5";
	gp+="set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$' offset 1";
	gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc="+tostring(bc_)+" $\\delta="+tostring(param)+"$'";
	gp+="plot '"+filename_+"-corr.dat' u 1:($6==1?$2:1/0):3 w errorbars lt 1 lc 1 lw 2 t 'Independant measures',\\";
	gp+="     '"+filename_+"-corr.dat' u 1:($6==0?$2:1/0):3 w errorbars lt 1 lc 2 lw 2 t 'Mean',\\";
	gp+="     "+tostring(poly_e(N_/m_-1)) + " w l lc 3 t 'd-merization="+tostring(poly_e(N_/m_-1)-poly_e(N_/m_-2))+"',\\";
	gp+="     "+tostring(poly_e(N_/m_-2)) + " w l lc 3 notitle";
	gp.save_file();
	rst.link_figure(analysis_+path_+dir_+filename_+"-corr.png","Correlation on links",analysis_+path_+dir_+filename_+"-corr.gp",1000);

	unsigned int length(long_range_corr_.size());
	if(length>0){
		Gnuplot gp(analysis_+path_+dir_,filename_+"-long-range-corr");
		gp+="stats '"+filename_+"-long-range-corr.dat' nooutput";
		gp.xrange(0,length+1);
		gp.yrange("1.1*STATS_min_y","1.1*STATS_max_y");
		gp+="set xlabel '$\\|i-j\\|$' offset 0,0.5";
		gp+="set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(j)>$' offset 1";
		gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$ bc="+tostring(bc_)+" $\\delta="+tostring(param)+"$'";
		gp+="set key right bottom";
		gp+="a=1.0";
		gp+="b=1.0";
		gp+="eta=1.0";
		gp+="m="+tostring(m_)+".0";
		gp+="N="+tostring(N_)+".0";
		gp+="f(x) = a/(x*x) + b*cos(2.0*pi*x*m/N)/(x**eta)";
		gp+="set fit quiet";
		switch(N_/m_){
			case 2:{ gp+="fit [3:"+tostring(length)+"] f(x) '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; } break;
			case 3:{
					   switch((length + 1) % 3){
						   case 0:{ gp+="fit [2:"+tostring(length)+"] f(x) '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
						   case 1:{ gp+="fit [5:"+tostring(length)+"] f(x) '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
						   case 2:{ gp+="fit [3:"+tostring(length)+"] f(x) '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
					   }break;
				   }break;
			default :{ gp+="fit ["+tostring(N_-1)+":] f(x) '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
		}
		gp+="plot for [IDX=0:"+tostring(nruns-1)+"] '"+filename_+"-long-range-corr.dat' i IDX u 1:($4==1?$2:1/0):3 w errorbars lt 1 lc 3 ps 0 notitle,\\";
		gp+="     for [IDX=0:"+tostring(nruns-1)+"] '"+filename_+"-long-range-corr.dat' i IDX u 1:($4==0?$2:1/0):3 w errorbars lt 1 lc 1 ps 0 notitle,\\";
		gp+="                   '"+filename_+"-long-range-corr.dat' i "+tostring(nruns)+" u 1:2:3 w errorbars lt 1 lc 2 lw 2 notitle,\\";
		gp+="                   f(x) notitle";
		gp.save_file();
		rst.link_figure(analysis_+path_+dir_+filename_+"-long-range-corr.png","Long range correlation",analysis_+path_+dir_+filename_+"-long-range-corr.gp",1000);
	}
	/*}*/

	rst.text(read.get_header());
	rst.save(false);
	all_link_names_.append(tostring(param));
}

void AnalyseParameter::extract_level_4(){
	/*comparison of E(param_optimal)|n=fixé*/
	IOFiles read(sim_+path_+dir_+filename_+".jdbin",false);

	unsigned int nof(0);
	double param(0.0);
	double polymerization_strength;
	read>>nof;

	Data<double> min_E;
	min_E.set_x(0.0);
	double min_param(0.0);
	double min_polymerization_strength(0.0);

	unsigned int idx(0);
	for(unsigned int i(0);i<nof;i++){
		read>>ref_>>N_>>m_>>n_>>bc_>>param>>E_>>polymerization_strength;
		if(E_.get_x()<min_E.get_x()){ 
			idx = i;
			min_E = E_;
			min_param = param;
			min_polymerization_strength = polymerization_strength;
		}
	}

	(*write_)("ref (type of wavefunction)",ref_);
	(*write_)("N (N of SU(N))",N_);
	(*write_)("m (# of particles per site)",m_);
	(*write_)("n (# of site)",n_);
	(*write_)("bc (boundary condition)",bc_);
	(*write_)("param",param);
	(*write_)("min E",min_E);
	(*write_)("polymerization strength",min_polymerization_strength);

	Gnuplot gp(analysis_+path_+dir_,filename_);
	gp+="set xlabel '$\\delta$' offset 0,1";
	gp+="set ylabel '$\\dfrac{E}{n}$' rotate by 0 offset 1";
	if(idx==0){
		gp+="f(x) = a+b*x**eta";
		gp+="a="+tostring(min_E.get_x());
		gp+="b=1";
		gp+="eta=1";
		gp+="set fit quiet";
		gp+="fit f(x) '"+filename_+".dat' u 1:($4==0?$2:1/0):3 via a,b,eta";
		gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$'";
		gp+="plot '"+filename_+".dat' u 1:($4==1?$2:1/0):3 w e t 'Independant measures',\\";
		gp+="     '"+filename_+".dat' u 1:($4==0?$2:1/0):3 w e t 'Mean',\\";
		gp+="     f(x) t sprintf('eta %3.4f',eta)";
	} else {
		gp+="f(x) = a+b*(x-c)*(x-c)";
		gp+="a="+tostring(min_E.get_x());
		gp+="b=1";
		gp+="c="+tostring(min_param);
		gp+="set fit quiet";
		gp+="fit f(x) '"+filename_+".dat' u 1:($4==0?$2:1/0):3 via a,b,c";
		gp+="set title '$N="+tostring(N_)+"$ $m="+tostring(m_)+"$ $n="+tostring(n_)+"$'";
		gp+="plot '"+filename_+".dat' u 1:($4==1?$2:1/0):3 w e t 'Independant measures',\\";
		gp+="     '"+filename_+".dat' u 1:($4==0?$2:1/0):3 w e t 'Mean',\\";
		gp+="     f(x) t sprintf('min %3.4f',c)";
	}
	gp.save_file();
	gp.create_image(true);
	all_link_names_.append(filename_);
}

void AnalyseParameter::extract_level_3(){
	/*evolution in function of n*/
	IOFiles read(sim_+path_+dir_+filename_+".jdbin",false);

	unsigned int nof;
	Vector<unsigned int> ref;
	unsigned int N;
	unsigned int m;
	unsigned int n;
	int bc;
	double param;
	Data<double> E;
	double polymerization_strength;
	read>>nof;
	for(unsigned int i(0);i<nof;i++){
		read>>ref>>N>>m>>n>>bc>>param>>E>>polymerization_strength;
		(*data_write_)<<n<<" "<<polymerization_strength<<IOFiles::endl;
		(*write_)("ref",ref);
		(*write_)("N",N);
		(*write_)("m",m);
		(*write_)("n",n);
		(*write_)("bc",bc);
		(*write_)("E",E);
		(*write_)("param",param);
		(*write_)("polymerization strength",polymerization_strength);
	}
	all_link_names_.append(filename_);
}

void AnalyseParameter::extract_level_2(){
	Gnuplot gp(analysis_+path_+dir_,filename_);
	gp+="set xlabel '$n^{-1}$'";
	gp.xrange(0,"");
	gp.yrange(0,"");
	gp+="set key bottom";
	gp+="plot '"+filename_+".dat' u (1/$1):2 t 'd-merization strength'";
	gp.save_file();
	gp.create_image(true);
	all_link_names_.append(filename_);
}
