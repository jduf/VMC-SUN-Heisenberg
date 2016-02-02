#ifndef DEF_KAGOMEDIRAC
#define DEF_KAGOMEDIRAC

#include "Kagome.hpp"

template<typename Type>
class KagomeDirac: public Kagome<Type>{
	public:
		KagomeDirac(System const& s);
		~KagomeDirac() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void display_results();

		Matrix<double> set_ab() const;
		unsigned int match_pos_in_ab(Vector<double> const& x) const;

		std::string extract_level_7();
		std::string extract_level_6();
		std::string extract_level_4();
};

template<typename Type>
KagomeDirac<Type>::KagomeDirac(System const& s):
	System(s),
	Kagome<Type>(set_ab(),6,"kagome-dirac")
{
	if(this->status_==2){
		this->init_fermionic();

		this->system_info_.text("KagomeDirac : 3 sites per unit cell, pi-flux per hexagon,");
		this->system_info_.text("no flux per triangle");
	}
}

/*{method needed for running*/
template<typename Type>
void KagomeDirac<Type>::compute_H(){
	//double t(1.0);
	this->H_.set(this->n_,this->n_,0);
	//Matrix<int> nb;
	//unsigned int s(0);
	//for(unsigned int i(0);i<this->Lx_;i++){
		//for(unsigned int j(0);j<this->Ly_;j++){
			///*site 0*/
			//s = this->spuc_*(i + j*this->Lx_);
			///*0-1*/this->H_(s,nb(0,0)) = nb(0,1)*t;
			///*0-1*/this->H_(s,nb(2,0)) = nb(2,1)*t;
//
			///*site 1*/
			//s++;
			///*0-1*/this->H_(s,nb(1,0)) = nb(1,1)*t;
			///*0-1*/this->H_(s,nb(3,0)) = -nb(3,1)*t;
//
			///*site 2*/
			//s++;
			///*0-1*/this->H_(s,nb(0,0)) = -nb(0,1)*t;
			///*0-1*/this->H_(s,nb(2,0)) = nb(2,1)*t;
//
			///*site 3*/
			//s++;
			///*0-1*/this->H_(s,nb(0,0)) = nb(0,1)*t;
			///*0-1*/this->H_(s,nb(2,0)) = -nb(2,1)*t;
//
			///*site 4*/
			//s++;
			///*0-1*/this->H_(s,nb(1,0)) = nb(1,1)*t;
			///*0-1*/this->H_(s,nb(3,0)) = -nb(3,1)*t;
//
			///*site 5*/
			//s++;
			///*0-1*/this->H_(s,nb(0,0)) = nb(0,1)*t;
			///*0-1*/this->H_(s,nb(2,0)) = nb(2,1)*t;
		//}
	//}
	this->H_ += this->H_.transpose();
}

template<typename Type>
Matrix<double> KagomeDirac<Type>::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2.0;
	tmp(1,0) = 0.0;
	tmp(0,1) =-1.0;
	tmp(1,1) = sqrt(3.0);
	return tmp;
}

template<typename Type>
unsigned int KagomeDirac<Type>::match_pos_in_ab(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match)){ return 0; }
	match(0) = 0.5;
	if(my::are_equal(x,match)){ return 1; }
	match(1) = 0.5;
	if(my::are_equal(x,match)){ return 2; }
	return 3;
}
/*}*/

/*{method needed for checking*/
template<typename Type>
void KagomeDirac<Type>::display_results(){
	//Matrix<int> nb;
	//double x0;
	//double x1;
	//double y0;
	//double y1;
	//double ll(1.0);
	//double ex(4.0*ll);
	//double exy(2.0*ll*cos(2.0*M_PI/6.0));
	//double ey(2.0*ll*sin(2.0*M_PI/6.0));
	//std::string color;
//
	//PSTricks ps("./","lattice");
	//ps.add("\\begin{pspicture}(-1,-1)(16,10)%"+this->filename_);
	//Matrix<double> cell(4,2);
	//cell(0,0) = 0.0;
	//cell(0,1) = 0.0;
	//cell(1,0) = ex;
	//cell(1,1) = 0.0;
	//cell(2,0) = ex + exy;
	//cell(2,1) = ey;
	//cell(3,0) = exy;
	//cell(3,1) = ey;
	//ps.polygon(cell,"linewidth=1pt,linecolor=red");
	//cell(1,0)*=this->Lx_;
	//cell(2,0) = this->Lx_*ex + this->Ly_*exy;
	//cell(2,1)*=this->Ly_;
	//cell(3,0)*=this->Ly_;
	//cell(3,1)*=this->Ly_;
	//ps.polygon(cell,"linewidth=1pt,linecolor=red,linestyle=dashed");

	//unsigned int s;
	//for(unsigned int i(0);i<this->Lx_;i++) {
		//for(unsigned int j(0);j<this->Ly_;j++) {
	//}
	//ps.add("\\end{pspicture}");
	//ps.save(true,true,true);
}

template<typename Type>
void KagomeDirac<Type>::check(){
}
/*}*/

/*{method needed for analysing*/
template<typename Type>
std::string KagomeDirac<Type>::extract_level_7(){
	this->rst_file_ = new RSTFile(this->info_+this->path_+this->dir_,this->filename_);
	unsigned int nruns;
	unsigned int tmax;

	(*this->read_)>>nruns>>tmax;
	this->obs_.push_back(Observable(*this->read_));
	this->jd_write_->write("E",this->obs_[0][0]);

	this->rst_file_->text(this->read_->get_header());
	this->rst_file_->save(false,true);
	delete this->rst_file_;
	this->rst_file_ = NULL;

	return this->filename_;
}

template<typename Type>
std::string KagomeDirac<Type>::extract_level_6(){
	unsigned int nof(0);
	(*this->read_)>>nof;
	this->save_input(*this->jd_write_);
	for(unsigned int i(0);i<nof;i++){
		(*this->read_)>>this->obs_[0][0];
		(*this->data_write_)<<this->M_(0)<<" "<<this->obs_[0][0]<<" "<<this->ref_(0)<<this->ref_(1)<<this->ref_(2)<<IOFiles::endl;
	}
	this->jd_write_->write("E",this->obs_[0][0]);

	return this->filename_;
}

template<typename Type>
std::string KagomeDirac<Type>::extract_level_4(){
	(*this->read_)>>this->obs_[0][0];
	this->jd_write_->write("E",this->obs_[0][0]);
	(*this->data_write_)<<this->M_(0)<<" "<<this->obs_[0][0]<<" "<<this->ref_(0)<<this->ref_(1)<<this->ref_(2)<<IOFiles::endl;

	return this->filename_;
}
/*}*/
#endif
