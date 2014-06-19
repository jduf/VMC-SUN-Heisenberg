#ifndef DEF_BANDSTRUCTURE
#define DEF_BANDSTRUCTURE

#include "Gnuplot.hpp"
#include "Lapack.hpp"

template<typename Type>
class BandStructure{
	public:
		BandStructure(Matrix<Type> T, Matrix<Type> P);
		BandStructure(Matrix<Type> T, Matrix<Type> Px, Matrix<Type> Py);

	protected:
		Matrix<Type> TP_;//!< translation operator along x-axis 
		Matrix<std::complex<double> > evec_;
		Vector<std::complex<double> > eval_;

		/*!Evaluate the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, unsigned int idx);
};
	
template<typename Type>
BandStructure<Type>::BandStructure(Matrix<Type> T, Matrix<Type> P):
	TP_(T+P)
{
	Lapack<Type> ES(&TP_,false,'G');
	ES.eigensystem(&eval_,&evec_);

	IOFiles w("spectrum.dat",true);
	for(unsigned int i(0);i<T.row();i++){
		w<<log(projection(P,i)).imag()<<" "<<projection(T,i).real()<<IOFiles::endl;
	}

	Gnuplot gp("./","spectrum");
	gp.xrange("-pi","pi");
	gp+="plot 'spectrum.dat' u 1:2";
	gp.save_file();
	gp.create_image(true);
}
	
template<typename Type>
BandStructure<Type>::BandStructure(Matrix<Type> T, Matrix<Type> Px, Matrix<Type> Py):
	TP_(T+Type(3.)*Px+Type(7.)*Py)
{
	//std::cout<<T*Px-Px*T<<std::endl;
	//std::cout<<T*Py-Py*T<<std::endl;
	
	Lapack<Type> ES(&TP_,false,'G');
	ES.eigensystem(&eval_,&evec_);
	IOFiles w("spectrum.dat",true);
	for(unsigned int i(0);i<T.row();i++){
		w<<log(projection(Px,i)).imag()<<" "<<log(projection(Py,i)).imag()<<" "<<projection(T,i).real()<<IOFiles::endl;
	}

	Gnuplot gp("./","spectrum");
	gp.xrange("-pi","pi");
	gp+="splot 'spectrum.dat' u 1:2:3";
	gp.save_file();
	gp.create_image(true);
}

template<typename Type>
std::complex<double> BandStructure<Type>::projection(Matrix<Type> const& O, unsigned int idx){
	Vector<std::complex<double> > tmp(O.row(),0.0);
	std::complex<double> out(0.0);
	for(unsigned int i(0);i<O.row();i++){
		for(unsigned int j(0);j<O.col();j++){
			tmp(i) += O(i,j)*evec_(j,idx);
		}
	}
	for(unsigned int i(0);i<O.row();i++){
		out += std::conj(evec_(i,idx))*tmp(i);
	}
	return out;
}
#endif
