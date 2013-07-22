#include "SquareFermi.hpp"

SquareFermi::SquareFermi(Parseur& P):
	Square<double>(P)
{
	if(!P.status()){
		if(N_site==N_row*N_col){
			compute_EVec();
			for(unsigned int spin(0);spin<N_spin;spin++){
				for(unsigned int i(0);i<N_site;i++){
					for(unsigned int j(0);j<N_m;j++){
						EVec(i+spin*N_site,j) = T(i,j);
					}
				}
			}
			if(successful){
				std::string filename("square-fermi");
				filename += "-N" + tostring(N_spin);
				filename += "-S" + tostring(N_site);
				filename += "-" + tostring(N_row) + "x" + tostring(N_col);
				if(bc == 1){ filename += "-P";} 
				else { filename += "-A";}
				save(filename);
			} else {
				std::cerr<<"SquareFermi : degeneate"<<std::endl;
			}
		} else {
			std::cerr<<"SquareFermi : the cluster is not a square"<<std::endl;
		}
	}
}

SquareFermi::~SquareFermi(){}

void SquareFermi::compute_EVec(){
	T.set(N_site,N_site,0.0);
	double t(-1.0);
	for(unsigned int i(0); i < N_site; i++){
		/*horizontal hopping*/
		if( (i+1) % N_col ){ T(i,i+1) = t;}	
		else { T(i+1-N_col,i) = bc*t;}
		/*vertical hopping*/
		if( i+N_col < N_site ){  T(i,i+N_col) = t; } 
		else { T(i-(N_row-1)*N_col,i) = bc*t;}
	}
	/*\warning if i take the transpose, the diagonal will be counted twice*/
	diagonalize_EVec('S');
}

void SquareFermi::save(std::string filename){
	Write w(filename+".jdbin");
	RST rst;
	rst.text("fermi : all colors experience the same Hamiltonian");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("is_complex",false);
	w("N_spin",N_spin);
	w("N_m",N_m);
	w("sts",sts);
	w("EVec",EVec);
	w("bc",bc);
	w("N_row",N_row);
	w("N_col",N_col);
}

