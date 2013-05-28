#include "CreateSystem.hpp"

template<>
void CreateSystem<double>::compute_T(unsigned int N_row, unsigned int N_col, double bc){
	H.set(0.0);
	T.set(0.0);
	is_complex = false;
	mat_type = 'S';
	double t(-1.0);
	switch(N_n){
		case 3: 
			{
				double th(1.0);
				double td(-1.0);
				//{//for SU(3) with 12 sites
					//unsigned int i(0);
					//unsigned int j(0);
					//for(unsigned int ix(0);ix<N_col;ix++){
						//i = j;
						//j = N_row*(2*ix+1)*6;
						//for(unsigned int iy(0);iy<N_row;iy++){ 
							//// cell (A)
							//H(i,j) = t;
							//H(i,i+1) = t;
							//H(j,j+1) = t;
							//T(i,j) = th;
							//T(i,i+1) = th;
							//T(j,j+1) = th;
							//i++;
							//j++;
							//if( j+N_row*6>N_site ){
								//H(j%(N_row*6),j)= bc*t; 
								//T(j%(N_row*6),j)= bc*td; 
							//} else { 
								//H(j,j+N_row*6) = t; 
								//T(j,j+N_row*6) = td; 
							//}
							//H(i,i+1) = t;
							//H(j,j+1) = t;
							//T(i,i+1) = th;
							//T(j,j+1) = th;
							//i++;
							//j++;			
							//// cell (B)
							//H(i,j) = t;
							//H(i,i+1) = t;
							//H(j,j+1) = t;
							//T(i,j) = th;
							//T(i,i+1) = td;
							//T(j,j+1) = td;
							//i++;
							//j++;
							//if( j+N_row*6>N_site ){
								//H(j%(N_row*6),j)= bc*t; 
								//T(j%(N_row*6),j)= bc*th; 
							//} else {
								//H(j,j+N_row*6) = t; 
								//T(j,j+N_row*6) = th; 
							//}
							//H(i,i+1) = t;
							//H(j,j+1) = t;
							//T(i,i+1) = th;
							//T(j,j+1) = th;
							//i++;
							//j++;			
							//// cell (C)
							//H(i,j) = t;
							//H(i,i+1) = t;
							//H(j,j+1) = t;
							//T(i,j) = td;
							//T(i,i+1) = th;
							//T(j,j+1) = th;
							//i++;
							//j++;
							//if( j+N_row*6>N_site ){
								//H(j%(N_row*6),j)= bc*t; 
								//T(j%(N_row*6),j)= bc*th; 
							//} else{
								//H(j,j+N_row*6) = t;
								//T(j,j+N_row*6) = th;
							//}
							//if( iy+1 == N_row  ){ // link (A)-(C)
								//H(i+1,j)= bc*t; 
								//H((i+1)-N_row*6,i) = bc*t; 
								//T(i+1,j)= bc*td; 
								//T((i+1)-N_row*6,i) = bc*td; 
							//} else { // link (A)-(C)
								//H(i,i+1) = t; 
								//H(j,j+1) = t; 
								//T(i,i+1) = td; 
								//T(j,j+1) = td; 
							//}
							//i++;
							//j++;
						//}
					//}
				//}
				//{//for SU(3) (0pipi) with 6 sites
					//unsigned int i(0);
					//for(unsigned int l(0);l<N_row;l++){
						//for(unsigned int c(0);c<N_col;c++){
							////std::cout<<i<<std::endl;
							//H(i,i+1) = t; //0-1
							//T(i,i+1) = th; //0-1
							//i++;//1
							//H(i,i+1) = t; //1-2
							//T(i,i+1) = th; //1-2
							//i++;//2
							//H(i,i+1) = t; //2-3
							//T(i,i+1) = td; //2-3
							//H(i,i+2) = t; //2-4
							//T(i,i+2) = th; //2-4
							//i++;//3
							//H(i,i+2) = t; //3-5
							//T(i,i+2) = th; //3-5
							//if(c+1==N_col){
								//H((i+3) - N_col*6, i) = t; // 0-3 b
								//T((i+3) - N_col*6, i) = bc*th; // 0-3 b
							//} else {
								//H(i,i+3) = t; //3-6
								//T(i,i+3) = th; //3-6
							//}
							//i++;//4
							//if(l+1==N_row){
								//H(i-4-l*N_col*6,i) = t; // 0-4 b
								//T(i-4-l*N_col*6,i) = bc*td; // 0-4 b
							//} else {
								//H(i,i-4+N_col*6) = t; // 4-18
								//T(i,i-4+N_col*6) = td; // 4-18
							//}
							//i++;//5
							//if(l+1==N_row){
								//H(i-4-l*N_col*6,i) = t; // 1-5 b
								//T(i-4-l*N_col*6,i) = bc*td; // 1-5 b
							//} else {
								//H(i,i-4+N_col*6) = t; // 5-19
								//T(i,i-4+N_col*6) = td; // 5-19
							//}
							////std::cout<<i<<std::endl;
							//if(c+1==N_col){
								//H((i+5) - N_col*6, i) = t; // 4-5 b
								//T((i+5) - N_col*6, i) = bc*th; // 4-5 b
							//} else {
								//H(i,i+5) = t; //5-10
								//T(i,i+5) = th; //5-10
							//}
							//i++;//6
						//}
					//}
				//}
				//{//for SU(3) (0pipi) with 6 sites
					//unsigned int i(0);
					//for(unsigned int l(0);l<N_row;l++){
						//for(unsigned int c(0);c<N_col;c++){
							////std::cout<<i<<std::endl;
							//H(i,i+1) = t; //0-1
							//T(i,i+1) = th; //0-1
							//i++;//1
							//H(i,i+1) = t; //1-2
							//T(i,i+1) = th; //1-2
							//i++;//2
							//H(i,i+1) = t; //2-3
							//T(i,i+1) = td; //2-3
							//H(i,i+2) = t; //2-4
							//T(i,i+2) = th; //2-4
							//i++;//3
							//H(i,i+2) = t; //3-5
							//T(i,i+2) = th; //3-5
							//if(c+1==N_col){
								//H((i+3) - N_col*6, i) = t; // 0-3 b
								//T((i+3) - N_col*6, i) = bc*th; // 0-3 b
							//} else {
								//H(i,i+3) = t; //3-6
								//T(i,i+3) = th; //3-6
							//}
							//i++;//4
							//if(l+1==N_row){
								//H(i-4-l*N_col*6,i) = t; // 0-4 b
								//T(i-4-l*N_col*6,i) = bc*td; // 0-4 b
							//} else {
								//H(i,i-4+N_col*6) = t; // 4-18
								//T(i,i-4+N_col*6) = td; // 4-18
							//}
							//i++;//5
							//if(l+1==N_row){
								//H(i-4-l*N_col*6,i) = t; // 1-5 b
								//T(i-4-l*N_col*6,i) = bc*td; // 1-5 b
							//} else {
								//H(i,i-4+N_col*6) = t; // 5-19
								//T(i,i-4+N_col*6) = td; // 5-19
							//}
							////std::cout<<i<<std::endl;
							//if(c+1==N_col){
								//H((i+5) - N_col*6, i) = t; // 4-5 b
								//T((i+5) - N_col*6, i) = bc*th; // 4-5 b
							//} else {
								//H(i,i+5) = t; //5-10
								//T(i,i+5) = th; //5-10
							//}
							//i++;//6
						//}
					//}
				//}
				{//for SU(4) (pipipi) with 4 sites
					unsigned int i(0);
					for(unsigned int l(0);l<N_row;l++){
						for(unsigned int c(0);c<N_col;c++){
							//std::cout<<i<<std::endl;
							H(i,i+1) = t; //0-1
							T(i,i+1) = th;
							H(i,i+3) = t; //0-3
							T(i,i+3) = th;
							i++;//1
							H(i,i+1) = t; //1-2
							T(i,i+1) = th; 
							i++;//2
							if(c+1==N_col){
								H(i+1-c*4, i) = t; // 2-3 b
								T(i+1-c*4, i) = bc*th;  
							} else {
								H(i,i+5) = t; 
								T(i,i+5) = th;
							}
							if(l+1==N_row){
								H(i-1-l*N_col*4,i) = t; // 2-1 b
								T(i-1-l*N_col*4,i) = bc*td; 
							} else {
								H(i,i-1+N_col*4) = t; 
								T(i,i-1+N_col*4) = td;
							}
							i++;//3
							if(l+1==N_row){
								H(i-3-l*N_col*4, i) = t; // 3-0 b
								T(i-3-l*N_col*4, i) = bc*th; 
							} else {
								H(i,i-3+N_col*4) = t; 
								T(i,i-3+N_col*4) = th;
							}
							i++;//4
						}
					}
				}
				//{//for SU(3) si cellule pi au milieu (vérifier que tout fonctionne) with 6 sites
					//unsigned int i(0);
					//for(unsigned int l(0);l<N_row;l++){
						//for(unsigned int c(0);c<N_col;c++){
							////std::cout<<i<<std::endl;
							////0
							//H(i,i+1) = t; //0-1
							//T(i,i+1) = td; 
							//i++;//1
							//H(i,i+1) = t; //1-2
							//T(i,i+1) = th; 
							//if(l != 0 && l+1 != N_row && c+1 != N_col){
								//H(i-9,i) = t; //1-4
								//T(i-9,i) = th; 
							//} else {
								//if(c+1==N_col){
									//if(l==0){
										//H(i,N_site-2-c*6) = t; // 1-4 bc bc
										//T(i,N_site-2-c*6) = th; // 1-4 bc bc
									//} else {
										//H(i+3-(2*N_col-1)*6,i) = t; //1-4
										//T(i+3-(2*N_col-1)*6,i) = th; 
									//}
								//} else {
									//if(l==0){
										//H(i,i+3+((N_row-1)*N_col+1)*6) = t; //1-4 bc
										//T(i,i+3+((N_row-1)*N_col+1)*6) = bc*th; 
									//} else {//l+1==N_row
										//H(i-9,i) = t; //1-4
										//T(i-9,i) = th; 
									//}
//
								//}
							//}
							//i++;//2
							////std::cout<<i<<std::endl;
							//H(i,i+1) = t; //2-3
							//T(i,i+1) = td;
							//if(c+1 == N_col){
								//H(i+3-c*6,i) = t; //2-5 bc
								//T(i+3-c*6,i) = bc*th;
							//} else {
								//H(i,i+9) = t;//2-5
								//T(i,i+9) = th;
							//}
							//i++;//3
							//H(i,i+1) = t; //3-4
							//T(i,i+1) = th;
							//if(l+1==N_row){
								//H(i-3-l*N_col*6,i) = t;//3-0 bc
								//T(i-3-l*N_col*6,i) = bc*th;
							//} else {
								//H(i,i-3+N_col*6) = t;//3-0
								//T(i,i-3+N_col*6) = th;
							//}
							//i++;//4
							//H(i,i+1) = t; //4-5
							//T(i,i+1) = td; 
							//i++;//5
							//H(i-5,i) = t; //0-5
							//T(i-5,i) = th; 
							////std::cout<<i<<std::endl;
							//i++;//6
						//}
					//}
				//}
				H += H.transpose();
				T += T.transpose();
				//std::cout<<T<<std::endl;
				//T.print_mathematica();
				break;
			}
		case 4:
			{
				//for(unsigned int i(0); i< N_row; i++){
					//for(unsigned int j(0); j< N_col; j++){
						//if(j+1 == N_col){ // x hopping
							//H(i*N_col , i*N_col + j) = t;
							//T(i*N_col , i*N_col + j) = bc*t;
						//} else {
							//H( i*N_col + j , i*N_col + j + 1) = t; 
							//T( i*N_col + j , i*N_col + j + 1) = t; 
						//}
						//if(i+1 == N_row ){// y hopping
							//H(j, i*N_col + j) = t;
							//T(j, i*N_col + j) = bc*t;
						//} else{
							//H(i*N_col + j, (i+1)*N_col + j) = t;
							//T(i*N_col + j, (i+1)*N_col + j) = t;
						//}
					//}
				//}
				//H += H.transpose();
				//T += T.transpose();
	// new csl : Vishwanath
				for(unsigned int i(0); i< N_row; i++){
					for(unsigned int j(0); j< N_col; j++){
						if(j+1 == N_col){// x hopping
							H(i*N_col , i*N_col + j) = t;
							if(i % 2 == 0){
								T(i*N_col , i*N_col + j) = bc*t;
							} else {
								T(i*N_col , i*N_col + j) = -bc*t;
							}
						} else {
							H( i*N_col + j , i*N_col + j + 1) = t; 
							if(i % 2 == 0){
								T( i*N_col + j , i*N_col + j + 1) = t; 
							} else {
								T( i*N_col + j , i*N_col + j + 1) = -t; 
							}
						}
						if(i+1 == N_row ){// y hopping
							H(j, i*N_col + j) = t;
							T(j, i*N_col + j) = bc*t;
						} else{
							H(i*N_col + j, (i+1)*N_col + j) = t;
							T(i*N_col + j, (i+1)*N_col + j) = t;
						}
					}
				}
				H += H.transpose();
				T += T.transpose();
				break;
			} 
		default:
			{
				std::cerr<<"CreateSystem : compute_T(N_row,N_col,bc) bad N_n"<<std::endl;
			}
	}
	//std::cout<<H<<std::endl;
}

template<>
void CreateSystem<std::complex<double> >::compute_T(unsigned int N_row, unsigned int N_col, double bc){
	std::cout<<N_site<<" "<<N_col<<" "<<N_row<<" "<<bc<<std::endl;
	H.set(0.0);
	T.set(0.0);
	is_complex = true;
	mat_type = 'H';
	double t(-1.0);
	double phi(2*M_PI/N_spin);
	unsigned int j(0);
	for(unsigned int i(0); i< N_row; i++){
		for(unsigned int j(0); j< N_col; j++){
			if(j+1 == N_col){
				H( i*N_col , i*N_col + j) = bc*t;
				T( i*N_col , i*N_col + j) = bc*t;
			} else { 
				H( i*N_col + j , i*N_col + j + 1) = t; 
				T( i*N_col + j , i*N_col + j + 1) = t; 
			}
			if(i+1 == N_row ){
				H(j, i*N_col + j) = bc*t;
				T(j, i*N_col + j) = bc*t*std::polar(1.0,-((j%N_spin)+1)*phi);
			} else{
				H(i*N_col + j, (i+1)*N_col + j) = t;
				T(i*N_col + j, (i+1)*N_col + j) = t*std::polar(1.0,((j%N_spin)+1)*phi);
			}
		}
	}
	H += H.transpose();
	T += T.trans_conj(); 
}

template<>
void CreateSystem<double>::compute_T(){
	is_complex = false;
	mat_type = 'S';
	double t(-1.0);
	for(unsigned int i(0); i< N_site-1; i++){
		H(i,i+1) = t;
		T(i,i+1) = t;
	}
	H += H.transpose();
	T += T.transpose();
}

template<>
void CreateSystem<std::complex<double> >::compute_T(){
	std::cerr<<"CreateSystem : compute_T not defined for a complex without specifing N_row and N_col"<<std::endl;
}
