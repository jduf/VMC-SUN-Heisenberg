#ifndef DEF_KAGOMEFERMI
#define DEF_KAGOMEFERMI

#include "Kagome.hpp"

template<typename Type>
class KagomeFermi: public Kagome<Type>{
	public:
		KagomeFermi(System const& s);
		~KagomeFermi() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void display_results();
		void lattice();

		Matrix<double> set_ab() const;
		unsigned int match_pos_in_ab(Vector<double> const& x) const;

		std::string extract_level_7();
		std::string extract_level_6();
		std::string extract_level_4();
};

template<typename Type>
KagomeFermi<Type>::KagomeFermi(System const& s):
	System(s),
	Kagome<Type>(set_ab(),3,"kagome-fermi")
{
	if(this->status_==2){
		this->init_fermionic();

		this->system_info_.text("KagomeFermi : All hopping term are identical, no flux, 3 sites per unit cell");
	}
}

/*{method needed for running*/
template<typename Type>
void KagomeFermi<Type>::compute_H(){
	this->H_.set(this->n_,this->n_,0);
	for(unsigned int i(0);i<this->obs_[0].nlinks(); i++){
		this->H_(this->obs_[0](i,0),this->obs_[0](i,1)) = this->obs_[0](i,4);
	}
	this->H_ += this->H_.transpose();
}

template<typename Type>
unsigned int KagomeFermi<Type>::match_pos_in_ab(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match)){ return 0; }
	match(0) = 0.5;
	if(my::are_equal(x,match)){ return 1; }
	match(1) = 0.5;
	if(my::are_equal(x,match)){ return 2; }
	return 3;
}

template<typename Type>
Matrix<double> KagomeFermi<Type>::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2.0;
	tmp(1,0) = 0.0;
	tmp(0,1) =-1.0;
	tmp(1,1) = sqrt(3.0);
	return tmp;
}
/*}*/

/*{method needed for checking*/
template<typename Type>
void KagomeFermi<Type>::lattice(){
	compute_H();

	Matrix<int> nb;
	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	Type t;
	PSTricks ps(this->info_+this->path_+this->dir_,this->filename_);
	ps.begin(-20,-10,20,10,this->filename_);

	Matrix<double> polygon(4,2);
	polygon(0,0)=0;
	polygon(0,1)=0;
	polygon(1,0)=this->LxLy_(0,0);
	polygon(1,1)=this->LxLy_(1,0);
	polygon(2,0)=this->LxLy_(0,0)+this->LxLy_(0,1);
	polygon(2,1)=this->LxLy_(1,0)+this->LxLy_(1,1);
	polygon(3,0)=this->LxLy_(0,1);
	polygon(3,1)=this->LxLy_(1,1);
	ps.polygon(polygon,"linecolor=green");

	polygon(0,0)=0;
	polygon(0,1)=0;
	polygon(1,0)=this->ab_(0,0);
	polygon(1,1)=this->ab_(1,0);
	polygon(2,0)=this->ab_(0,0)+this->ab_(0,1);
	polygon(2,1)=this->ab_(1,0)+this->ab_(1,1);
	polygon(3,0)=this->ab_(0,1);
	polygon(3,1)=this->ab_(1,1);
	ps.polygon(polygon,"linecolor=black");

	for(unsigned int i(0);i<this->n_;i++){
		xy0 = this->get_pos_in_lattice(i);
		this->set_pos_LxLy(xy0);
		xy0 = (this->LxLy_*xy0).chop();
		ps.put(xy0(0)-0.20,xy0(1)+0.15,my::tostring(i));

		nb = this->get_neighbourg(i);

		t = this->H_(i,nb(0,0));
		if(std::abs(t)>1e-4){
			xy1 = this->get_pos_in_lattice(nb(0,0));
			this->set_pos_LxLy(xy1);
			xy1 = this->LxLy_*xy1;
			if((xy0-xy1).norm_squared()>1.0001){
				linestyle = "dashed"; 
				xy1 = xy0;
				if(i%3==2){
					xy1(0) += this->dir_nn_(2,0);
					xy1(1) += this->dir_nn_(2,1);
				} else {
					xy1(0) += this->dir_nn_(0,0);
					xy1(1) += this->dir_nn_(0,1);
				}
				ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(0,0)));
			} else { linestyle = "solid";  }
			xy1 = xy1.chop();

			if(my::real(t)>0){ color = "blue";}
			else { color = "red"; }
			linewidth = my::tostring(std::abs(t))+"mm";
			/*(+x)-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}

		t = this->H_(i,nb(1,0));
		if(std::abs(t)>1e-4){
			xy1 = this->get_pos_in_lattice(nb(1,0));
			this->set_pos_LxLy(xy1);
			xy1 = this->LxLy_*xy1;
			if((xy0-xy1).norm_squared()>1.0001){
				linestyle = "dashed"; 
				xy1 = xy0;
				xy1(0) += this->dir_nn_(1,0);
				xy1(1) += this->dir_nn_(1,1);
				ps.put(xy1(0)-0.20,xy1(1)+0.15,my::tostring(nb(1,0)));
			} else { linestyle = "solid";  }
			xy1 = xy1.chop();

			if(my::real(t)>0){ color = "blue";}
			else { color = "red"; }
			linewidth = my::tostring(std::abs(t))+"mm";
			/*(+y)-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}

	}
	ps.end(true,true,true);
}

template<typename Type>
void KagomeFermi<Type>::check(){
	///*{debug 1*/
	//Matrix<int> nb;
	//for(unsigned int i(0);i<this->n_;i++){
	//nb = this->get_neighbourg(i);
	//std::cout<<i<<" ";
	//for(unsigned int j(0);j<z_;j++){
	//std::cout<<nb(j,0)<<" ";
	//}
	//std::cout<<std::endl;
	//}
	///*}*/
	///*{debug 2*/
	//Matrix<int> nb;
	//double t(1.0);
	//Matrix<double> Ttest(this->n_,this->n_,0);
	//for(unsigned int s(0);s<this->n_;s++){
	//nb = this->get_neighbourg(s);
	//for(unsigned int i(0);i<z_;i++){ Ttest(s,nb(i,0)) = t; }
	//}
	//for(unsigned int i(0);i<this->n_;i++){
	//for(unsigned int j(0);j<this->n_;j++){
	//if(std::abs(Ttest(i,j)-this->H_(i,j))>0.2){
	//std::cout<<i<<" "<<j<<std::endl;
	//}
	//}
	//}
	///*}*/
	///*{debug 3*/
	//unsigned int k(0);
	//for(unsigned int i(0);i<this->n_;i++){
	//for(unsigned int j(0);j<this->n_;j++){
	//if(this->H_(i,j)!=0){
	//k++;
	////std::cout<<i<<" "<<j<<" "<<this->H_(i,j)<<std::endl;
	//}
	//}
	//}
	//std::cout<<k<<" "<<links_.row()<<std::endl;
	///*}*/
	///*{debug 4*/
	//Matrix<int> nb;
	//for(unsigned int s(0);s<this->n_;s++){
	//nb = this->get_neighbourg(s);
	//for(unsigned int i(0);i<z_;i++){
	//if(nb(i,1)<0){ std::cout<<s<<" "<<nb(i,0)<<std::endl; }
	//}
	//}
	///*}*/

	//this->plot_band_structure();

	///*{Rotation matrix : [R,T]!=0*/
	//Matrix<Type> R(this->H_.row(),this->H_.col(),0);
	//Matrix<int> nb;
	//unsigned int s;
	//unsigned int rs;
	//int bc;
	//for(unsigned int j(0); j<this->Ly_;j++){
	//for(unsigned int i(0); i<this->Lx_;i++){
	//bc = 1;
	//s = this->spuc_*(i + j*this->Lx_);
	//
	//nb = this->get_neighbourg(s);
	//rs = nb(0,0);
	//bc*= nb(0,1);
	//nb = this->get_neighbourg(rs);
	//rs = nb(1,0);
	//bc*= nb(1,1);
	//for(unsigned int k(0); k<2*i;k++){
	//nb = this->get_neighbourg(rs);
	//rs = nb(1,0);
	//bc*= nb(1,1);
	//}
	//if(j!=0){
	//nb = this->get_neighbourg(rs);
	//rs = nb(1,0);
	//bc*= nb(1,1);
	//for(unsigned int k(0); k<2*j-1;k++){
	//nb = this->get_neighbourg(rs);
	//rs = nb(2,0);
	//bc*= nb(2,1);
	//}
	//nb = this->get_neighbourg(rs);
	//rs = nb(3,0);
	//bc*= nb(3,1);
	//}
	//nb = this->get_neighbourg(rs);
	//R(s,rs) = bc;
	//R(s+1,nb(0,0)) = bc*nb(0,1);
	//R(s+2,nb(1,0)) = bc*nb(1,1);
	////std::cout<<s<<" "<<rs<<" "<<bc<<std::endl;
	////std::cout<<s+1<<" "<<nb(0,0)<<" "<<bc*nb(0,1)<<std::endl;
	////std::cout<<s+2<<" "<<nb(1,0)<<" "<<bc*nb(1,1)<<std::endl;
	//}
	//}
	////std::cout<<R*R*R*R*R*R<<std::endl;
	////std::cout<<R*this->H_-this->H_*R<<std::endl;
	////std::cout<<R*R*R*this->Tx_-this->Tx_*R*R*R<<std::endl;
	///*}*/
	//this->plot_band_structure();
	this->info_ = "";
	this->path_ = "";
	this->dir_  = "./";
	this->filename_ ="kagome-fermi";
	display_results();
}

template<typename Type>
void KagomeFermi<Type>::display_results(){
	lattice();
}
/*}*/

/*{method needed for analysing*/
template<typename Type>
std::string KagomeFermi<Type>::extract_level_7(){
	this->rst_file_ = new RSTFile(this->info_+this->path_+this->dir_,this->filename_);
	unsigned int nruns;
	unsigned int tmax;

	(*this->read_)>>nruns>>tmax;
	this->save_input(*this->jd_write_);
	this->save_output(*this->jd_write_);

	this->rst_file_->text(this->read_->get_header());
	this->rst_file_->save(false,true);
	delete this->rst_file_;
	this->rst_file_ = NULL;

	return this->filename_;
}

template<typename Type>
std::string KagomeFermi<Type>::extract_level_6(){
	unsigned int nof(0);
	(*this->read_)>>nof;
	for(unsigned int i(0);i<nof;i++){
		(*this->data_write_)<<this->M_(0)<<" "<<this->E_.get_x()<<" "<<this->E_.get_dx()<<" "<<this->ref_(0)<<this->ref_(1)<<this->ref_(2)<<IOFiles::endl;
	}
	this->save_input(*this->jd_write_);
	this->save_output(*this->jd_write_);

	return this->filename_;
}

template<typename Type>
std::string KagomeFermi<Type>::extract_level_4(){
	(*this->read_)>>this->E_;
	this->jd_write_->write("energy per site",this->E_);
	(*this->data_write_)<<this->M_(0)<<" "<<this->E_.get_x()<<" "<<this->E_.get_dx()<<" "<<this->ref_(0)<<this->ref_(1)<<this->ref_(2)<<IOFiles::endl;

	return this->filename_;
}
/*}*/
#endif
