#ifndef DEF_KAGOME
#define DEF_KAGOME

#include "System2D.hpp"

template<typename Type>
class Kagome: public System2D<Type>{
	public:
		/*{Description*/
		/*!Constructor that organises the n sites according to the ratio Lx/Ly
		 * for a system with spuc sites per unit cell. Calls the
		 * GenericSystem<Type>(4,filename), to construct a system with 4 links
		 * per sites */
		/*}*/
		Kagome(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Kagome()=0;

	protected:
		void set_observables(unsigned int const& which);
		Vector<double> get_pos_in_lattice(unsigned int const& i) const { return Vector<double>(i); }
		unsigned int match_pos_in_ab(Vector<double> const& x) const { (void)(x); return 0; }

	private:
		Matrix<double> set_geometry(unsigned int const& n) const;
		Vector<double> vector_towards(unsigned int const& i, unsigned int const& dir) const { return Vector<double>(i,dir); }
		void try_neighbourg(Vector<double>& tn, unsigned int const& j) const { (void)(tn); (void)(j); }
};

template<typename Type>
Kagome<Type>::Kagome(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(set_geometry(this->n_),Matrix<double>(Lx,Ly),spuc,4,filename)
{
	std::cerr<<__PRETTY_FUNCTION__<<" : new def of set_nn_links will be problematic"<<std::endl;
	if(this->status_==2){ this->set_nn_links(Vector<unsigned int>(1,2)); }
}

template<typename Type>
Kagome<Type>::~Kagome() = default;

template<typename Type>
void Kagome<Type>::set_observables(unsigned int const& which){
	this->E_.set(50,5,false);
	this->corr_types_.resize(which);

	if(which>0){
		this->corr_types_[0].set(this->link_types_[0].row(),50,5,false);
	}
	if(which>1){ /*the long range correlation*/
		this->corr_types_[1].set(this->n_,50,5,false);
		this->link_types_.push_back(Matrix<int>(this->n_,2));
		for(unsigned int i(0);i<this->n_;i++){
			this->link_types_[1](i,0) = 0;
			this->link_types_[1](i,1) = i;
		}
	}
	if(which==6){ /*the (anti)symmetric correlation*/
		this->corr_types_[2].set(this->n_/2,50,5,false);
		this->corr_types_[3].set(this->n_/2,50,5,false);
		this->corr_types_[4].set(this->n_/2,50,5,false);
		this->corr_types_[5].set(this->n_/2,50,5,false);

		this->link_types_.push_back(Matrix<int>(this->n_/2,2));
		this->link_types_.push_back(Matrix<int>(this->n_/2,2));
		this->link_types_.push_back(Matrix<int>(this->n_/2,2));
		this->link_types_.push_back(Matrix<int>(this->n_/2,2));
		for(unsigned int i(0);i<this->n_/2;i++){
			this->link_types_[2](i,0) = 0;
			this->link_types_[2](i,1) = 2*i;
			this->link_types_[3](i,0) = 0;
			this->link_types_[3](i,1) = 2*i+1;
			this->link_types_[4](i,0) = 1;
			this->link_types_[4](i,1) = 2*i;
			this->link_types_[5](i,0) = 1;
			this->link_types_[5](i,1) = 2*i+1;
		}
	}
}

//template<typename Type>
//Matrix<int> Kagome<Type>::get_neighbourg(unsigned int i) const {
//Matrix<int> nb(this->z_,2,1);
//switch(this->spuc_){
//case 3: case 6:
//{
//switch(i%3){
//case 0:
//{
///*+x neighbour*/
//nb(0,0) = i+1;
///*+x+y neighbour*/
//nb(1,0) = i+2;
///*-x neighbour*/
//if(i%(this->spuc_*this->Lx_)){ nb(2,0) = i-2; }
//else {
//nb(2,0) = i-2+this->spuc_*this->Lx_;
//nb(2,1) = this->bc_;
//}
///*-x-y neighbour*/
//if(i+1>=this->spuc_*this->Lx_){ nb(3,0) = i+2-this->spuc_*this->Lx_; } 
//else {
//nb(3,0) = i+2+(this->Ly_-1)*this->spuc_*this->Lx_;
//nb(3,1) = this->bc_;
//}
//}break;
//case 1:
//{
///*+x neighbour*/
//if((i+2)%(this->spuc_*this->Lx_)){ nb(0,0) = i+2; }
//else {
//nb(0,0) = i+2-this->spuc_*this->Lx_;
//nb(0,1) = this->bc_;
//}
///*-x+y neighbour*/
//nb(1,0) = i+1;
///*-x neighbour*/
//nb(2,0) = i-1;
///*+x-y neighbour*/
//if(i+1>=this->spuc_*this->Lx_){
//if((i+2)%(this->spuc_*this->Lx_)){ nb(3,0) = i+4-this->spuc_*this->Lx_; } 
//else {
//nb(3,0) = i+4-2*this->spuc_*this->Lx_; 
//nb(3,1) = this->bc_;
//}
//} else {
//if((i+2)%(this->spuc_*this->Lx_)){
//nb(3,0) = i+4+(this->Ly_-1)*this->spuc_*this->Lx_;
//nb(3,1) = this->bc_;
//} else {
//nb(3,0) = this->Lx_*(this->Ly_-1)*this->spuc_+2; 
//nb(3,1) = this->bc_*this->bc_;
//}
//}
//}break;
//case 2:
//{
///*+x+y neighbour*/
//if(i<this->n_-this->spuc_*this->Lx_){ nb(0,0) = i+this->spuc_*this->Lx_-2; }
//else { 
//nb(0,0) = i-2-this->spuc_*this->Lx_*(this->Ly_-1);
//nb(0,1) = this->bc_;
//}
///*-x+y neighbour*/
//if((i-2)%(this->spuc_*this->Lx_)){
//if(i<this->n_-this->spuc_*this->Lx_){ nb(1,0) = i-4+this->spuc_*this->Lx_; }
//else {
//nb(1,0) = i-4-this->spuc_*this->Lx_*(this->Ly_-1);
//nb(1,1) = this->bc_;
//}
//} else {
//if(i<this->n_-this->spuc_*this->Lx_){
//nb(1,0) = i-4+2*this->spuc_*this->Lx_;
//nb(1,1) = this->bc_;
//} else {
//nb(1,0) = this->spuc_*this->Lx_-2;
//nb(1,1) = this->bc_*this->bc_;
//}
//}
///*-x-y neighbour*/
//nb(2,0) = i-2;
///*+x-y neighbour*/
//nb(3,0) = i-1;
//}break;
//}
//}break;
//case 9:
//{
//switch(i%9){
//case 0:
//{ 
//nb(0,0) = i+1; 
//if(i%(this->Lx_*this->spuc_)){ nb(1,0) = i-1; }
//else { 
//nb(1,0) = i-1+this->Lx_*this->spuc_; 
//nb(1,1) = this->bc_; 
//}
//if(i>=(this->Lx_*this->spuc_)){ nb(2,0) = i+6-this->Lx_*this->spuc_; }
//else { 
//nb(2,0) = i+6+(this->Ly_-1)*this->Lx_*this->spuc_;
//nb(2,1) = this->bc_; 
//}
//nb(3,0) = i+5;
//}break;
//case 1:
//{
//nb(0,0) = i+1;
//if((i-1)%(this->Lx_*this->spuc_)!=0){ nb(1,0) = i-3; }
//else {
//nb(1,0) = i-3+this->Lx_*this->spuc_; 
//nb(1,1) = this->bc_; 
//}
//if((i-1)%(this->Lx_*this->spuc_)){ nb(2,0) = i-2; }
//else {
//nb(2,0) = i-2+this->Lx_*this->spuc_; 
//nb(2,1) = this->bc_; 
//}
//nb(3,0) = i-1;
//}break;
//case 2:
//{
//nb(0,0) = i+1;
//nb(1,0) = i+4;
//if((i-2)%(this->Lx_*this->spuc_)!=0){ nb(2,0) = i-4; }
//else {
//nb(2,0) = i-4+this->Lx_*this->spuc_; 
//nb(2,1) = this->bc_; 
//}
//nb(3,0) = i-1;
//}break;
//case 3:
//{
//nb(0,0) = i+1;
//nb(1,0) = i+5;
//nb(2,0) = i+3;
//nb(3,0) = i-1;
//}break;
//case 4:
//{
//nb(0,0) = i+1;
//if(i>=(this->Lx_*this->spuc_)){ nb(1,0) = i+3-this->Lx_*this->spuc_; }
//else { 
//nb(1,0) = i+3+(this->Ly_-1)*this->Lx_*this->spuc_;
//nb(1,1) = this->bc_; 
//}
//nb(2,0) = i+4;
//nb(3,0) = i-1;
//}break;
//case 5:
//{
//nb(0,0) = i-5;
//if(i>=(this->Lx_*this->spuc_)){ nb(1,0) = i+1-this->Lx_*this->spuc_; }
//else { 
//nb(1,0) = i+1+(this->Ly_-1)*this->Lx_*this->spuc_;
//nb(1,1) = this->bc_; 
//}
//if(i>=(this->Lx_*this->spuc_)){ nb(2,0) = i+2-this->Lx_*this->spuc_; }
//else { 
//nb(2,0) = i+2+(this->Ly_-1)*this->Lx_*this->spuc_;
//nb(2,1) = this->bc_; 
//}
//nb(3,0) = i-1;
//}break;
//case 6:
//{
//if(i<(this->Ly_-1)*this->Lx_*this->spuc_){ nb(0,0) = i-1+this->Lx_*this->spuc_; }
//else {
//nb(0,0) = i-1-(this->Ly_-1)*this->Lx_*this->spuc_;
//nb(0,1) = this->bc_; 
//}
//if(i<(this->Ly_-1)*this->Lx_*this->spuc_){ nb(1,0) = i-6+this->Lx_*this->spuc_; }
//else {
//nb(1,0) = i-6-(this->Ly_-1)*this->Lx_*this->spuc_;
//nb(1,1) = this->bc_; 
//}
//nb(2,0) = i-4;
//nb(3,0) = i-3;
//}break;
//case 7:
//{
//if((i+2)%(this->Lx_*this->spuc_)!=0){ nb(0,0) = i+3; }
//else {
//nb(0,0) = i+3-this->Lx_*this->spuc_;
//nb(0,1) = this->bc_; 
//}
//if((i+2)%(this->Lx_*this->spuc_)!=0){ nb(1,0) = i+4; }
//else {
//nb(1,0) = i+4-this->Lx_*this->spuc_;
//nb(1,1) = this->bc_; 
//}
//if(i<(this->Ly_-1)*this->Lx_*this->spuc_){ nb(2,0) = i-3+this->Lx_*this->spuc_; }
//else {
//nb(2,0) = i-3-(this->Ly_-1)*this->Lx_*this->spuc_;
//nb(2,1) = this->bc_; 
//}
//if(i<(this->Ly_-1)*this->Lx_*this->spuc_){ nb(3,0) = i-2+this->Lx_*this->spuc_; }
//else {
//nb(3,0) = i-2-(this->Ly_-1)*this->Lx_*this->spuc_;
//nb(3,1) = this->bc_; 
//}
//}break;
//case 8:
//{
//if((i+1)%(this->Lx_*this->spuc_)!=0){ nb(0,0) = i+1; }
//else {
//nb(0,0) = i+1-this->Lx_*this->spuc_;
//nb(0,1) = this->bc_; 
//}
//if((i+1)%(this->Lx_*this->spuc_)!=0){ nb(1,0) = i+2; }
//else {
//nb(1,0) = i+2-this->Lx_*this->spuc_;
//nb(1,1) = this->bc_; 
//}
//nb(2,0) = i-5;
//nb(3,0) = i-4;
//}break;
//}
//}break;
//}
//return nb;
//}

template<typename Type>
Matrix<double> Kagome<Type>::set_geometry(unsigned int const& n) const {
	Matrix<double> tmp(n,n);
	std::cerr<<__PRETTY_FUNCTION__<<" : undefined"<<std::endl;
	return tmp;
}
#endif
