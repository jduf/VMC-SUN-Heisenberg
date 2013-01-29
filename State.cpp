#include "State.hpp"

/*Constructors and destructor*/
/*{*/
State::State(System *S,bool whole_copy):
	S(S),
	A(new Matrice[S->N_spin]),
	Ainv(new Matrice[S->N_spin]),
	s(new unsigned int[S->N_site]),
	wis(new unsigned int[S->N_site]),
	det(0.0),
	delete_matrix(true),
	whole_copy(whole_copy)
{
	srand(time(NULL));
	unsigned int site(0);
	unsigned int N_as(S->N_site);
	unsigned int available_sites[S->N_site];
	for(unsigned int i(0); i < N_as; i++){
		available_sites[i]  = i;	
	}

	for(unsigned int i(0); i<S->N_spin; i++){
		for(unsigned int j(0); j<S->N_m; j++){
			site = rand() % N_as;
			s[i*S->N_m+j] = available_sites[site];
			wis[available_sites[site]] = i*S->N_m+j; 	
			for(unsigned int k(site); k < N_as-1; k++){
				available_sites[k] = available_sites[k+1];
			}
			N_as--;
		}
	}
	init_matrices(S->N_m,S->N_spin);
	//for(unsigned int i(0); i < S->N_spin; i++){
		//Lapack A_(A[i].ptr(),A[i].size(),'G');
		//A_.inv();
	//}
	compute_det();
	//std::cout<<"creation"<<std::endl;
}

State::State(State const& s):
	S(s.S),
	A(s.A),
	Ainv(s.Ainv),
	s(s.s),
	wis(s.wis),
	det(s.det),
	delete_matrix(false),
	whole_copy(false)
{
	//std::cout<<"copie"<<std::endl;
}

State::~State(){
	//std::cout<<"destructeur : state"<<std::endl;
	if(delete_matrix){
		delete[] A;
		delete[] Ainv;
	}
	delete[] s;
	delete[] wis;
}
/*}*/

/*operator*/
/*{*/
State& State::operator=(State const& s){
	//this->S = s.S; // inutile car pointe de toute façon sur le bon system
	
	if(whole_copy){
		for(unsigned int i(0); i < S->N_site; i++){ //attention, S->N_site doit != 0
			this->s[i] = s.s[i];
			this->wis[i] = s.wis[i];
		}
		for(unsigned int i(0); i< S->N_spin; i++){
			this->A[i] = s.A[i];
		}
	} else {
		for(unsigned int i(0);i<2;i++){
			this->column_changed[i] = s.column_changed[i];
			this->matrix_changed[i] = s.matrix_changed[i];
		}
	}
	//std::cout<<"affectation : state"<<std::endl;
	return (*this);
}

std::ostream& operator<<(std::ostream& flux, State const& S){
	S.color(flux);
	return flux;
}

double operator/(State const& Snew, State const& Sold){
	double d(1.0),w(0.0);
	unsigned int N(4);
	for(unsigned int i(0); i<2;i++){
		w=compute_ratio_det(mold,mnew,Snew.get_changed_column(i),N);
		std::cout<<lnew.det()/lold.det()<<" "<<w<<std::endl;
		d *= w;
	}
	return d;
}
/*}*/

/*methods that return something related to the class*/
/*{*/
State State::swap() const {
	unsigned int s1,s2,p1,p2,a,b;
	p1 = rand() % S->N_m;
	p2 = rand() % S->N_m;
	s1 = rand() % S->N_spin;
	s2 = rand() % S->N_spin;
	while(s1==s2){
		s2 = rand() % S->N_spin;
	}
    a = s[s1*S->N_m + p1];
	b = s[s2*S->N_m + p2];
	
	return swap(a,b);
}

State State::swap(unsigned int a, unsigned int b) const {
	State new_s(*this);
	
	new_s.s[wis[b]] = s[wis[a]];
	new_s.s[wis[a]] = s[wis[b]];

	new_s.wis[b] = wis[a];
	new_s.wis[a] = wis[b];

	new_s.modify_matrix(wis[a],wis[b]);

	return new_s;
}

void State::color(std::ostream& flux) const{
	unsigned int col(0);
	for(unsigned int i(0); i<S->N_site;i++){
		col = wis[i] / S->N_m;
		flux << col<<" ";
	}
}

unsigned int State::get_changed_column(unsigned int i) const {
	return column_changed[i];
}

double compute_ratio_det(Matrice const& mold, Matrice const& mnew, unsigned int col, unsigned int N){
	double d(0.0);
	for(unsigned int i(0);i<N;i++){
		d += mold(i,col)*mnew(i,col);
	}
	return d;
}
/*}*/

/*methods that modify the class*/
/*{*/
void State::init_matrices(unsigned int N_m, unsigned int N_spin){
	unsigned int l(0);
	for(unsigned int i(0); i < N_spin; i++){
		A[i] = Matrice (N_m);
		for(unsigned int j(0); j < N_m; j++){
			l = s[i*S->N_m+j];
			for(unsigned int k(0); k < N_m; k++){
				A[i](k,j) = S->U(l,k);
			}
		}
	}
}

void State::compute_det(){
	det = 1.0;
	for(unsigned int i(0); i<S->N_spin; i++){
		Lapack Ai(A[i],'G');
		det *= Ai.det();
	}

}


void State::modify_matrix(unsigned int wisa, unsigned int wisb){
	matrix_changed[0] = wisa / S->N_m;
	matrix_changed[1] = wisb / S->N_m;
	column_changed[0] = wisa % S->N_m;
	column_changed[1] = wisb % S->N_m;
}
/*}*/

/*other methods*/
/*{*/
void State::print() const{
	//std::cout<<"{";
	//for(unsigned int i(0); i<S->N_spin; i++){
		//std::cout<<"{ ";
		//for(unsigned int j(0); j<S->N_m; j++){
			//std::cout<<s[i*S->N_m+j]<<" ";
		//}
		//std::cout<<"}";
	//}
	//std::cout<<"}"<<std::endl;
	//for(unsigned int i(0); i<S->N_spin; i++){
		//for(unsigned int j(0); j<S->N_m; j++){
			//std::cout<<wis[i*S->N_m+j]<<" ";
		//}
	//}
	//std::cout<<std::endl;
	for(unsigned int i(0);i<S->N_spin;i++){
		std::cout<<std::endl;
		A[i].print();
	}
	std::cout<<det<<std::endl<<std::endl;
}
/*}*/
