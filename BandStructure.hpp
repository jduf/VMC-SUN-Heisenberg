

template<typename Type>
class BandStructure{
	public:
		BandStructure(Matrix<Type> T, Matrix<Type> P);
		BandStructure(Matrix<Type> T, Matrix<Type> Px, Matrix<Type> Py);
		~Square();

	protected:
		Matrix<Type> TP_;//!< translation operator along x-axis 

		/*!Evaluate the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<Type> const& O, Matrix<std::complex<double> > const& base, unsigned int bra, unsigned int ket);

		void save(Vector<double> kx, Vector<double> ky, Vector<double> E);
		void save(Vector<double> k, Vector<double> E);
};
		
template<typename Type>
BandStructure<Type>::BandStructure(Matrix<Type> T, Matrix<Type> Px, Matrix<Type> Py):
	TP_(T+3.*Px+7*Py),
{
	//std::cout<<T*Px-Px*T<<std::endl;
	//std::cout<<T*Py-Py*T<<std::endl;
	//
	unsigned int n(T.row());

	Vector<std::complex<double> > eval;
	Matrix<std::complex<double> > evec;
	Lapack<double> ES(&TP_,false,'G');
	ES.eigensystem(&eval,&evec);
	Vector<double> kx(n_,1);
	Vector<double> ky(n_,1);
	Vector<double> E(n_,1);
	for(unsigned int i(0);i<n;i++){
		kx(i) = log(projection(Px_,evec,i,i)).imag();
		ky(i) = log(projection(Py_,evec,i,i)).imag();
		E(i) = projection(T,evec,i,i).real();
	}
	save(kx,ky,E);
}

template<typename Type>
BandStructure<Type>::BandStructure(Matrix<Type> T, Matrix<Type> P):
	TP_(T+P),
{
	//std::cout<<T*P-P*T<<std::endl;
	unsigned int n(T.row());

	Vector<std::complex<double> > eval;
	Matrix<std::complex<double> > evec;
	Lapack<double> ES(&TP,false,'G');
	ES.eigensystem(&eval,&evec);
	Vector<double> k(n_,1);
	Vector<double> E(n_,1);
	for(unsigned int i(0);i<n;i++){
		k(i) = log(projection(P,evec,i,i)).imag();
		E(i) = projection(T,evec,i,i).real();
	}
	save(k,E);
}

template<typename Type>
std::complex<double> BandStructure<Type>::projection(Matrix<Type> const& O, Matrix<std::complex<double> > const& base, unsigned int bra, unsigned int ket){
	Vector<std::complex<double> > tmp(O.row(),0.0);
	std::complex<double> out(0.0);
	for(unsigned int i(0);i<O.row();i++){
		for(unsigned int j(0);j<O.col();j++){
			tmp(i) += O(i,j)*base(j,ket);
		}
	}
	for(unsigned int i(0);i<O.row();i++){
		out += std::conj(base(i,bra))*tmp(i);
	}
	return out;
}

template<typename Type>
void BandStructure<Type>::save(Vector<double> kx, Vector<double> ky, Vector<double> E){
	Gnuplot gp(filename_+"-band-structure","splot");
	gp.save_data(filename_+"-spectrum",kx,ky,E);
	gp.add_plot_param(" ,\\\n");
	Vector<unsigned int> index(E.sort());
	gp.save_data(filename_+"-spectrum-sorted",kx.sort(index).range(0,m_),ky.sort(index).range(0,m_),E.range(0,m_));
}

template<typename Type>
void BandStructure<Type>::save(Vector<double> k, Vector<double> E){
	Gnuplot gp(filename_+"-band-structure","plot");
	gp.save_data(filename_+"-spectrum",kE);
	gp.add_plot_param(" ,\\\n");
	Vector<unsigned int> index(E.sort());
	gp.save_data(filename_+"-spectrum-sorted",k.sort(index).range(0,m_),E.range(0,m_));
}
