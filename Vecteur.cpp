#include "Vecteur.hpp"

/*Constructors and destructor*/
/*{*/
Vecteur::Vecteur():v(NULL),N(0){}

Vecteur::Vecteur(unsigned int N):
	v(new double[N]),
	N(N)
{
	std::cout<<"taille"<<std::endl;
	fill_vecteur(0.0);
}

Vecteur::Vecteur(unsigned int N, double val):
	v(new double[N]),
	N(N)
{
	std::cout<<"taille+const"<<std::endl;
	fill_vecteur(val);
}

Vecteur::Vecteur(Vecteur const& vec):
	v(new double[vec.size()]),
	N(vec.size())
{
	std::cout<<"copie"<<std::endl;
	for(unsigned int i(0);i<N;i++){
			v[i] = vec(i);
	}
}

Vecteur::~Vecteur(){
	delete v;
	std::cout<<"destructeur"<<std::endl;
}
/*}*/

/*operators*/
/*{*/
Vecteur& Vecteur::operator=(Vecteur const& vec){
	std::cout<<"affectation : vecteur "<<std::endl;
	if(this->N!=vec.N){
		std::cerr<<"impossible d'affecter deux Vecteurs si dim1 =! dim2"<<std::endl;
	}
	for(unsigned int i(0); i<N; i++){
		this->v[i] = vec.v[i];
	}
	return (*this);
}

double& Vecteur::operator()(unsigned int i){
	if(i >= N ){std::cout<<"bug taille"<<std::endl;}
	return v[i];
}

double const& Vecteur::operator()(unsigned int i) const{
	if(i >= N ){std::cout<<"bug taille"<<std::endl;}
	return v[i];
}

//Vecteur& Vecteur::operator*=(double const& d){
	//for(unsigned int i(0);i<N;i++){
		//v[i] *= d;
	//}
	//return (*this);
//}

//double operator*(Vecteur const& vec1, Vecteur const& vec2){
	//double s(0.0);
	//for(unsigned int i(0);i<N;i++){
		//s += vec1(i)*vec2(i);
	//}
	//return s;
//}
/*}*/

/*methods that modify the class*/
/*{*/
void  Vecteur::fill_vecteur(double val){
	for(unsigned int i(0);i<N;i++){
		v[i]=val;
	}
}

void Vecteur::chop(){
	for(unsigned int i(0);i<N;i++){
		if(std::fabs(v[i]) < 1e-10 ){v[i]=0;}
	}
}
/*}*/

/*other methods*/
/*{*/
void Vecteur::print() const{
	for(unsigned int i(0);i<N;i++){
		std::cout<<v[i]<<std::endl;;
	}
}

void Vecteur::test() const{
	print();
}
/*}*/
