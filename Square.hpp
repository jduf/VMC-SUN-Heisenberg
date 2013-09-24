#ifndef DEF_SQUARE
#define DEF_SQUARE

#include "CreateSystem.hpp"

template<typename Type>
class Square: public CreateSystem<Type>{
	public:
		Square(Parseur& P);
		~Square();

	protected:
		unsigned int Ly_, Lx_;
		Matrix<Type> Px_;
		Matrix<Type> Py_;

		void compute_H();
		void save(std::string filename);
		void compute_band_structure();
};

template<typename Type>
Square<Type>::Square(Parseur& P):
	CreateSystem<Type>(P,4),
	Ly_(floor(sqrt(this->n_))),
	Lx_(floor(sqrt(this->n_))),
	Px_(this->n_,this->n_,0.0),
	Py_(this->n_,this->n_,0.0)
{
	P.set("bc",this->bc_);
	if(!P.status()){
		if(this->n_==Ly_*Lx_){
			compute_H();
			this->compute_sts();
		} else {
			std::cerr<<"Square : the cluster is not a square"<<std::endl;
		}
	}
}

template<typename Type>
Square<Type>::~Square(){}

template<typename Type>
void Square<Type>::compute_H(){
	for(unsigned int i(0); i < this->n_; i++){
		/*horizontal hopping*/
		if( (i+1) % Lx_ ){ this->H_(i,i+1) = 1;}	
		else { this->H_(i+1-Lx_,i) = 1;}
		/*vertical hopping*/
		if( i+Lx_ < this->n_ ){  this->H_(i,i+Lx_) = 1; } 
		else { this->H_(i-(Ly_-1)*Lx_,i) = 1;}
	}
	this->H_ += this->H_.transpose();
}

template<typename Type>
void Square<Type>::compute_band_structure(){
	Matrix<double> TP(this->T_+3.*Px_+7.*Py_);
	Vector<std::complex<double> > eval;
	Matrix<std::complex<double> > evec;
	Lapack<double> ES(&TP,false,'G');
	ES.eigensystem(&eval,&evec);
	Vector<double> kx(this->n_,1);
	Vector<double> ky(this->n_,1);
	Vector<double> E(this->n_,1);
	for(unsigned int i(0);i<this->n_;i++){
		kx(i) = log(projection(Px_,evec,i,i)).imag();
		ky(i) = log(projection(Py_,evec,i,i)).imag();
		E(i) = projection(this->T_,evec,i,i).real();
	}
	Gnuplot gp("spectrum","3D");
	gp.save_data("spectrum",kx,ky,E);
	gp.code(" ,\\\n");
	Vector<unsigned int> index(E.sort());
	gp.save_data("spectrum-sorted",kx.sort(index).range(0,this->m_),ky.sort(index).range(0,this->m_),E.range(0,this->m_));
	gp.save_code();
}
#endif

