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
		bool is_over(){std::cout<<"bla"<<std::endl; return this->over_; }

	protected:
		void compute_H();
		void display_results();
		void lattice();

		std::string extract_level_7();
		std::string extract_level_6();
		std::string extract_level_4();
};

template<typename Type>
KagomeFermi<Type>::KagomeFermi(System const& s):
	System(s),
	Kagome<Type>(1,1,3,"kagome-fermi")
{
	if(this->status_==2){
		this->init_fermionic();

		this->system_info_.text("KagomeFermi : All hopping term are identical, no flux, 3 sites per unit cell");
	}
}

/*{method needed for running*/
template<typename Type>
void KagomeFermi<Type>::compute_H(){
	//double t(1.0);
	this->H_.set(this->n_,this->n_,0);
	//Matrix<int> nb;
	//unsigned int s(0);
	//for(unsigned int i(0);i<this->Lx_;i++){
		//for(unsigned int j(0);j<this->Ly_;j++){
			///*site 0*/
			//s = this->spuc_*(i + j*this->Lx_);
			//nb = this->get_neighbourg(s);
			///*0-1*/this->H_(s,nb(0,0)) = nb(0,1)*t;
			///*0-1*/this->H_(s,nb(2,0)) = nb(2,1)*t;
//
			///*site 1*/
			//s++;
			//nb = this->get_neighbourg(s);
			///*0-1*/this->H_(s,nb(1,0)) = nb(1,1)*t;
			///*0-1*/this->H_(s,nb(3,0)) = nb(3,1)*t;
//
			///*site 2*/
			//s++;
			//nb = this->get_neighbourg(s);
			///*0-1*/this->H_(s,nb(0,0)) = nb(0,1)*t;
			///*0-1*/this->H_(s,nb(2,0)) = nb(2,1)*t;
		//}
	//}
	this->H_ += this->H_.transpose();
}
/*}*/

/*{method needed for checking*/
template<typename Type>
void KagomeFermi<Type>::lattice(){
	//Matrix<int> nb;
	//double x0;
	//double x1;
	//double y0;
	//double y1;
	//double ll(1.0);
	//double ex(2.0*ll);
	//double exy(2.0*ll*cos(2.0*M_PI/6.0));
	//double ey(2.0*ll*sin(2.0*M_PI/6.0));
	//std::string color("black");
//
	//PSTricks ps(info_+path_+dir_,filename_+"-pstricks");
	//ps.begin(-1,-1,16,10,this->filename_);
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
//
	//unsigned int s;
	//for(unsigned int i(0);i<this->Lx_;i++) {
		//for(unsigned int j(0);j<this->Ly_;j++) {
			///*site 0*/
			//s = this->spuc_*(i+j*this->Lx_);
			//nb = this->get_neighbourg(s);
			//x0 = 0.2+i*ex+j*exy;
			///*0.05 is there so there is no problem with latex and it shows
			// * better which sites are in the unit cell*/
			//y0 = 0.1+j*ey;
			//ps.put(x0-0.2,y0+0.2,my::tostring(s));
			//x1 = x0+ll;
			//y1 = y0;
			//if(my::real(this->H_(s,nb(0,0)))>0){ color = "green"; }
			//else { color = "blue"; }
			///*0-1*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			//x1 = x0-ll;
			//if(my::real(this->H_(s,nb(2,0)))>0){ color = "green"; }
			//else { color = "blue"; }
			///*0-1*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
//
			///*site 1*/
			//s++;
			//nb = this->get_neighbourg(s);
			//x0 = x0+ll;
			//ps.put(x0+0.2,y0+0.2,my::tostring(s));
			//x1 = x0+ll*cos(4.0*M_PI/6.0);
			//y1 = y0+ll*sin(4.0*M_PI/6.0);
			//if(my::real(this->H_(s,nb(1,0)))>0){ color = "green"; }
			//else { color = "blue"; }
			///*1-2*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			//x1 = x0+ll*cos(10.0*M_PI/6.0);
			//y1 = y0+ll*sin(10.0*M_PI/6.0);
			//if(my::real(this->H_(s,nb(3,0)))>0){ color = "green"; }
			//else { color = "blue"; }
			///*1-2*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
//
			///*site 2*/
			//s++;
			//nb = this->get_neighbourg(s);
			//x0 = x0+ll*cos(4.0*M_PI/6.0);
			//y0 = y0+ll*sin(4.0*M_PI/6.0);
			//ps.put(x0+0.2,y0,my::tostring(s));
			//x1 = x0+ll*cos(2.0*M_PI/6.0);
			//y1 = y0+ll*sin(2.0*M_PI/6.0);
			//if(my::real(this->H_(s,nb(0,0)))>0){ color = "green"; }
			//else { color = "blue"; }
			///*2-0*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
			//x1 = x0+ll*cos(8.0*M_PI/6.0);
			//y1 = y0+ll*sin(8.0*M_PI/6.0);
			//if(my::real(this->H_(s,nb(2,0)))>0){ color = "green"; }
			//else { color = "blue"; }
			///*2-0*/	ps.line("->",x0,y0,x1,y1,"linewidth=1pt,linecolor="+color);
		//}
	//}
//
	//ps.end(true,true,true);
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

	this->plot_band_structure();

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
