#include "Honeycomb.hpp"

Honeycomb::Honeycomb(Parseur& P):
	CreateSystem<double>(P,3),
	N_row(floor(sqrt(N_m))),
	N_col(floor(sqrt(N_m)))
{
	P.set("bc",bc);
	if(!P.status()){
	if(N_m==N_row*N_col){
		mat_type='S';
		compute_T();
		compute_sts();
		compute_EVec();
		if(successful){
			std::string filename("honeycomb");
			filename +="-N" + tostring(N_spin);
			filename +="-S" + tostring(N_site);
			filename += "-" + tostring(N_row) +"x"+ tostring(N_col);
			if(bc == 1){ filename += "-P";} 
			else { filename += "-A";}
			save(filename);
		} else {
			if(bc == 1){
				std::cerr<<"CreateSystem : degeneate for PBC"<<std::endl;
			} else {
				std::cerr<<"CreateSystem : degeneate for APBC"<<std::endl;
			}
		}
	} else {
		std::cerr<<"CreateSystem : the cluster is not a square"<<std::endl;
	}
	}
}

Honeycomb::~Honeycomb(){}

void Honeycomb::compute_T(){
	double t(-1.0);
	double th(1.0);
	double td(t);
	//{//for SU(3) with 12 sites
	//unsigned int i(0);
	//unsigned int j(0);
	//for(unsigned int ix(0);ix<N_col;ix++){
	//i = j;
	//j = N_row*(2*ix+1)*6;
	//for(unsigned int iy(0);iy<N_row;iy++){ 
	////cell (A)
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
	////cell (B)
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
	////cell (C)
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
	//std::cout<<i<<std::endl;
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
	//std::cout<<i<<std::endl;
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
	{//for SU(4) (pipipi) with 4 sites
		if(N_spin!=4){std::cerr<<"works only with SU(4)"<<std::endl;}
		unsigned int i(0);
		for(unsigned int l(0);l<N_row;l++){
			for(unsigned int c(0);c<N_col;c++){
				H(i,i+1) = -1; //0-1
				T(i,i+1) = th;
				H(i,i+3) = -1; //0-3
				T(i,i+3) = th;
				i++;//1
				H(i,i+1) = -1; //1-2
				T(i,i+1) = th; 
				i++;//2
				if(c+1==N_col){
					H(i+1-c*4, i) = -1; // 2-3 b
					T(i+1-c*4, i) = bc*th;  
				} else {
					H(i,i+5) = -1; 
					T(i,i+5) = th;
				}
				if(l+1==N_row){
					H(i-1-l*N_col*4,i) = -1; // 2-1 b
					T(i-1-l*N_col*4,i) = bc*td; 
				} else {
					H(i,i-1+N_col*4) = -1; 
					T(i,i-1+N_col*4) = td;
				}
				i++;//3
				if(l+1==N_row){
					H(i-3-l*N_col*4, i) = -1; // 3-0 b
					T(i-3-l*N_col*4, i) = bc*th; 
				} else {
					H(i,i-3+N_col*4) = -1; 
					T(i,i-3+N_col*4) = th;
				}
				i++;//4
			}
		}
	}
	H += H.transpose();
	T += T.transpose();
}

void Honeycomb::save(std::string filename){
	Write w(filename+".jdbin");
	std::cerr<<"detail what kind of honeycomb it is"<<std::endl;
	RST rst;
	rst.text("Honeycomb ");
	rst.np();
	rst.title("Input values","~");

	w.set_header(rst.get());
	w("is_complex",false);
	w("N_spin",N_spin);
	w("N_m",N_m);
	w("sts",sts);
	w("EVec",T);
	w("bc",bc);
	w("N_row",N_row);
	w("N_col",N_col);
}

//{
//filename="honeycomb"+filename;
//Read r(hopfile);
//Matrice<double> tmp(N_spin*N_m);
//r>>tmp;
//double t(-1.0);
//double th(-1.0);
//for(unsigned int i(0);i<N_site;i++){
//for(unsigned int j(0);j<N_site;j++){
//if(std::abs(tmp(i,j)+1.0) < 1e-4){//original file th=-1
//H(i,j) = t;
//T(i,j) = th;
//}
//if(std::abs(tmp(i,j)-2.0) < 1e-4){ //original file td=2
//H(i,j) = t;
//T(i,j) = td;
//}
//if(std::abs(tmp(i,j)-1.0) < 1e-4){
//H(i,j) = t;
//T(i,j) = bc*th;
//}
//if(std::abs(tmp(i,j)+2.0) < 1e-4){
//H(i,j) = t;
//T(i,j) = bc*td;
//}
//}
//}
//mat_type = 'S';
//compute_sts();
//compute_EVec();
//if(successful){
//save();
//(*w)("th",th);
//(*w)("td",td);
//(*w)("bc",bc);
//} else {
//std::cerr<<"CreateSystem : degeneate"<<std::endl;
//}
//}
