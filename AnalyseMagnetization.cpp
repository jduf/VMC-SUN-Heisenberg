#include "AnalyseMagnetization.hpp"

void AnalyseMagnetization::open_files(std::string const& jdfile, std::string const& datafile, Directory const& d){
	if(level_>1){ 
		jd_write_ = new IOFiles(jdfile,true);
		(*jd_write_)("number of jdfiles",d.size());
		jd_write_->add_to_header("\n");
	}
	if(level_==5 || level_ == 4){
		data_write_ = new IOFiles(datafile,true);
	}
}

void AnalyseMagnetization::close_files(){
	if(jd_write_){ 
		switch(level_){
			case 5:{ rst_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","E.png",analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000); } break;
			case 4:{ rst_.last().link_figure(analyse_+path_+dir_.substr(0,dir_.size()-1)+".png","M.png",analyse_+path_+dir_.substr(0,dir_.size()-1)+".gp",1000); } break;
		}
		rst_.last().text(jd_write_->get_header());
		delete jd_write_;
		jd_write_ = NULL;
	}
	if(data_write_){
		delete data_write_;
		data_write_ = NULL;
	}
}

void AnalyseMagnetization::extract_level_5(){
	/*E(param)|n=fixé*/
	std::cerr<<"here should call theblalb"<<std::endl;
	//IOFiles read_(sim_+path_+dir_+filename_+".jdbin",false);
	//RSTFile rst(info_+path_+dir_,filename_);
//
	//unsigned int nruns;
//
	//read>>nruns>>ref_>>N_>>m_>>n_>>M_>>bc_;
	//data_write_->precision(10);
	//(*data_write_)<<"% M_0 E dE 0|1"<<IOFiles::endl;
	///* the +1 is the averages over all runs */
	//for(unsigned int i(0);i<nruns+1;i++){ 
		//read>>E_>>corr_>>long_range_corr_;
		//(*data_write_)<<M_(0)<<" "<<E_.get_x()<<" "<<E_.get_dx()<<" "<<(i<nruns?true:false)<<IOFiles::endl;
	//}
//
	//(*write_)("ref (type of wavefunction)",ref_);
	//(*write_)("N (N of SU(N))",N_);
	//(*write_)("m (# of particles per site)",m_);
	//(*write_)("n (# of site)",n_);
	//(*write_)("M (# of particles for each color)",M_);
	//(*write_)("bc (boundary condition)",bc_);
	//(*write_)("E",E_);
//
	//rst.text(read.get_header());
	//rst.save(false);
	//all_link_names_.append(tostring(M_(0)));
}

void AnalyseMagnetization::extract_level_4(){
	/*comparison of E(param_optimal)|n=fixé*/
	//IOFiles read(sim_+path_+dir_+filename_+".jdbin",false);
//
	//unsigned int nof(0);
	//read>>nof;
//
	//Vector<double> E(nof);
	//Vector<unsigned int> M(nof);
	//for(unsigned int i(0);i<nof;i++){
		//read>>ref_>>N_>>m_>>n_>>M_>>bc_>>E_;
		//E(i) = E_.get_x();
		//M(i) = M_(0);
	//}
	//Vector<unsigned int> index;
	//M.sort(std::greater<unsigned int>(),index);
	//E = E.order(index);
//
	//Vector<double> h(nof);
	//h(0) = 0;
	//for(unsigned int i(1);i<nof;i++){
		//h(i) = (E(i)-E(i-1))/(M(i-1)-M(i));
	//}
	//for(unsigned int i(0);i<nof;i++){
		//(*data_write_)<<h(i)<<" "<<double(M(0)-M(i))/M(0)<<IOFiles::endl;
	//}
//
	//Gnuplot gp(analyse_+path_+dir_,filename_);
	//gp+="plot '"+filename_+".dat' u 1:($4==1?$2:1/0):3 w e t 'Independant measures',\\";
	//gp+="     '"+filename_+".dat' u 1:($4==0?$2:1/0):3 w e t 'mean',\\";
	//gp.save_file();
	//gp.create_image(true);
	//all_link_names_.append(filename_);
}

void AnalyseMagnetization::extract_level_3(){
	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp+="plot '"+filename_+".dat' u 1:2 t 'Independant measures',\\";
	gp.save_file();
	gp.create_image(true);
	all_link_names_.append(filename_);
}
