#include "SquareMu.hpp"

SquareMu::SquareMu(Parseur& P):
	Square<double>(P),
	mu(P.get<double>("mu"))
{
	if(!P.status()){
		if(N_site==N_row*N_col){
			for(unsigned int spin(0);spin<N_spin;spin++){
				compute_EVec(spin);
				for(unsigned int i(0);i<N_site;i++){
					for(unsigned int j(0);j<N_m;j++){
						EVec(i+spin*N_site,j) = T(i,j);
					}
				}
			}
			if(successful){
				std::string filename("square-stripe");
				filename += "-N" + tostring(N_spin);
				filename += "-S" + tostring(N_site);
				filename += "-" + tostring(N_row) + "x" + tostring(N_col);
				if(bc == 1){ filename += "-P";} 
				else { filename += "-A";}
				filename += "-mu" + tostring(mu);
				save(filename);
			} else {
				std::cerr<<"SquareMu : degeneate"<<std::endl;
			}
		} else {
			std::cerr<<"SquareMu : the cluster is not a square"<<std::endl;
		}
	}
}

SquareMu::~SquareMu(){}

void SquareMu::compute_EVec(unsigned int spin){
	T.set(N_site,N_site,0.0);
	double t(-1.0);
	for(unsigned int i(0); i < N_site; i++){
		/*chemical potential*/
		if( (i-spin) % N_spin == 0 && i >= spin){ T(i,i) = mu; }
		/*horizontal hopping*/
		if( (i+1) % N_col ){ T(i,i+1) = t;}	
		else { T(i+1-N_col,i) = bc*t;  spin++; }
		/*vertical hopping*/
		if( i+N_col < N_site ){  T(i,i+N_col) = t; } 
		else { T(i-(N_row-1)*N_col,i) = bc*t;}
	}
	/*\warning if I take the transpose, the diagonal will be counted twice*/
	//T += T.transpose();
	//show(T,spin%N_spin+1);
	diagonalize_EVec('S');
}

void SquareMu::save(std::string filename){
	Write w(filename+".jdbin");
	RST rst;
	rst.text("Stripe order : each color lives on its own sublattice");
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
	w("mu",mu);
}

void SquareMu::show(Matrix<double> const& T,unsigned int spin){
	std::cout<<"T="<<std::endl;
	std::cout<<T<<std::endl;
	std::cout<<"favored sites :"<<std::endl;
	for(unsigned int i(0);i<N_row;i++){
		for(unsigned int j(0);j<N_col;j++){
			if(T(i+j*N_row,i+j*N_row)!=0){
				std::cout<<spin<<" ";
			} else { 
				std::cout<<0<<" ";
			}
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}
