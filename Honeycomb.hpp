#ifndef DEF_HONEYCOMB
#define DEF_HONEYCOMB

#include "CreateSystem.hpp"

template<typename Type>
class Honeycomb: public CreateSystem<Type>{
	public:
		Honeycomb(Parseur& P, std::string filename);
		virtual ~Honeycomb();

	protected:
		unsigned int Lx_;//!< dimension of the lattice along x-axis
		unsigned int Ly_;//!< dimension of the lattice along y-axis
		Matrix<Type> Px_;//!< translation operator along x-axis 
		Matrix<Type> Py_;//!< translation operator along y-axis 

		virtual Vector<unsigned int> get_neighbourg(unsigned int i);
};

template<typename Type>
Honeycomb<Type>::Honeycomb(Parseur& P, std::string filename):
	CreateSystem<Type>(P,3,filename),
	Lx_(floor(sqrt(this->m_))),
	Ly_(floor(sqrt(this->m_)))
{
	this->ref_(0) = 6;
	std::cerr<<"HoneycombSU3 will rewrite all the honeycomb things as a rectangular 12 sites unit cell"<<std::endl;
	this->bc_= P.get<double>("bc");
	if(!P.status()){
		if(this->m_==Ly_*Lx_){
			this->compute_sts();
			this->filename_ +="-N" + tostring(this->N_);
			this->filename_ +="-S" + tostring(this->n_);
			this->filename_ += "-" + tostring(Lx_) +"x"+ tostring(Ly_);
		} else {
			std::cerr<<"Honeycomb<Type> : the cluster is not a square"<<std::endl;
		}
	}
}

template<typename Type>
Honeycomb<Type>::~Honeycomb(){
	std::cerr<<"attention n!=Nm => System won't work"<<std::endl;
}

template<typename Type>
Vector<unsigned int> Honeycomb<Type>::get_neighbourg(unsigned int i){
	std::cerr<<"honecomb : need to define get_neighbourg "<<i<<std::endl;

	//unsigned int i(0);
	//for(unsigned int l(0);l<Ly_;l++){
	//for(unsigned int c(0);c<Lx_;c++){
	////0
	//this->H_(i,i+1) = 1;
	//if(l+1<Ly_){
	//this->H_(i,i+1+Lx_*4) = 1;
	//} else {
	//this->H_(i+1-l*Lx_*4,i) = 1;
	//}
	//if(c==0){
	//this->H_(i,i+Lx_*4-1) = 1;
	//} else {
	//this->H_(i-1,i) = 1;
	//}
	//i+=2;//2
	//this->H_(i,i-1) = 1;
	//this->H_(i,i+1) = 1; 
	//if(l==0){
	//this->H_(i,i+1+(Ly_-1)*Lx_*4) = 1;
	//} else {
	//this->H_(i,i+1-Lx_*4) = 1;
	//}
	//i+=2;//4
	//}
	//}
	//this->H_ += this->H_.transpose();
	return 0;
}
#endif
