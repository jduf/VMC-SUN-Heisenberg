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
	//std::cout<<"creation"<<std::endl;
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
	compute_det();
	for(unsigned int i(0);i<2;i++){
		cc[i] = 0;
		mc[i] = 0;
	}
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
	//std::cout<<"destructeur : state";
	if(delete_matrix){
	//std::cout<<" effectif";
		delete[] A;
		delete[] Ainv;
		delete[] s;
		delete[] wis;
	}
	//std::cout<<std::endl;
}
/*}*/

/*operator*/
/*{*/
State& State::operator=(State const& s){
	//this->S = s.S; // inutile car pointe de toute façon sur le bon system
	for(unsigned int i(0);i<2;i++){
		this->cc[i] = s.cc[i];
		this->mc[i] = s.mc[i];
	}

	if(whole_copy){
		//std::cout<<"whole ";
		unsigned int a(0),b(0);
		a = s.mc[0]*S->N_m+s.cc[0];
		b = s.mc[1]*S->N_m+s.cc[1];
		std::cout<<a<<" ab "<<b<<std::endl;
		std::cout<<s.s[a]<<" sab "<<s.s[b]<<std::endl;
		std::cout<<this->s[a]<<" sab "<<this->s[b]<<std::endl;
		
		this->s[b] = s.s[a];
		this->s[a] = s.s[b];

		std::cout<<this->s[a]<<" sab "<<this->s[b]<<std::endl;
		a = s.s[a];
		b = s.s[b];

		this->wis[b] = s.wis[a];
		this->wis[a] = s.wis[b];

		for(unsigned int i(0); i<S->N_m; i++){
			this->A[s.mc[0]](i,s.cc[0]) = s.A[s.mc[1]](i,s.cc[1]);
			this->A[s.mc[1]](i,s.cc[1]) = s.A[s.mc[0]](i,s.cc[0]);
		}

		for(unsigned int i(0); i<S->N_spin;i++){
			this->Ainv[i] = this->A[i];
			Lapack A_(this->Ainv[i].ptr(),this->Ainv[i].size(),'G');
			A_.inv();
		}
		this->compute_det();
	} 
	//std::cout<<"affectation"<<std::endl;
	return (*this);
}

std::ostream& operator<<(std::ostream& flux, State const& S){
	S.color(flux);
	return flux;
}

double operator/(State const& Snew, State const& Sold){
	return Snew.divided();
}
/*}*/

/*methods that return something related to the class*/
/*{*/
State State::swap() const {
	unsigned int s1,s2;
	s1 = rand() % S->N_spin;
	s2 = rand() % S->N_spin;
	while(s1==s2){
		s2 = rand() % S->N_spin;
	}
	
	State new_s(*this);

	new_s.mc[0] = s1;
	new_s.mc[1] = s2;
	new_s.cc[0] = rand() % S->N_m;
	new_s.cc[1] = rand() % S->N_m;
	
	//std::cout<<"swap"<<std::endl;
	//std::cout<<s[s1*S->N_m+new_s.cc[0]]<<" "<<s[s2*S->N_m+new_s.cc[1]]<<std::endl;

	return new_s;
}

State State::swap(unsigned int a, unsigned int b) const {
	State new_s(*this);

	new_s.mc[0] = wis[a] / S->N_m;
	new_s.mc[1] = wis[b] / S->N_m;
	new_s.cc[0] = wis[a] % S->N_m;
	new_s.cc[1] = wis[b] % S->N_m;

	return new_s;
}

void State::color(std::ostream& flux) const{
	unsigned int col(0);
	for(unsigned int i(0); i<S->N_site;i++){
		col = wis[i] / S->N_m;
		flux << col<<" ";
	}
}

double State::divided() const {
	double d1(0.0),d2(0.0);
	for(unsigned int i(0);i<S->N_m;i++){
		d1 += Ainv[mc[0]](cc[0],i)*A[mc[1]](i,cc[1]);
		d2 += Ainv[mc[1]](cc[1],i)*A[mc[0]](i,cc[0]);
	}
	return d1*d2;
}
/*}*/

/*methods that modify the class*/
/*{*/
void State::init_matrices(unsigned int N_m, unsigned int N_spin){
	unsigned int l(0);
	for(unsigned int i(0); i < N_spin; i++){
		A[i] = Matrice (N_m);
		Ainv[i] = Matrice (N_m);
		for(unsigned int j(0); j < N_m; j++){
			l = s[i*S->N_m+j];
			for(unsigned int k(0); k < N_m; k++){
				A[i](k,j) = S->U(l,k);
			}
		}
		Ainv[i] = A[i];
		Lapack A_(Ainv[i].ptr(),Ainv[i].size(),'G');
		A_.inv();
	}
}

void State::compute_det(){
	det = 1.0;
	for(unsigned int i(0); i<S->N_spin; i++){
		Lapack Ai(A[i],'G');
		det *= Ai.det();
	}

}
/*}*/

/*other methods*/
/*{*/
void State::print() const{
	std::cout<<"{";
	for(unsigned int i(0); i<S->N_spin; i++){
		std::cout<<"{ ";
		for(unsigned int j(0); j<S->N_m; j++){
			std::cout<<s[i*S->N_m+j]<<" ";
		}
		std::cout<<"}";
	}
	std::cout<<"}"<<std::endl;
	for(unsigned int i(0); i<S->N_site; i++){
		std::cout<<wis[i]<<" ";
	}
	std::cout<<std::endl;
	for(unsigned int i(0);i<2;i++){
		std::cout<<mc[i]<<" "<< cc[i]<<std::endl;
	}
	for(unsigned int i(0);i<S->N_spin;i++){
		std::cout<<std::endl;
		A[i].print();
	}
	for(unsigned int i(0);i<S->N_spin;i++){
		std::cout<<std::endl;
		(A[i]*Ainv[i]).print();
	}
}
/*}*/
