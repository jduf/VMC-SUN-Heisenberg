#include "State.hpp"

/*Constructors and destructor*/
/*{*/
State::State(System *S):
	S(S),
	A(new Matrice[S->N_spin]),
	s(new unsigned int[S->N_site]),
	wis(new unsigned int[S->N_site]),
	det(0.0)
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
	init_A(S->N_m,S->N_spin);
	compute_matrices();
	compute_det();
	//std::cout<<"creation"<<std::endl;
}

State::State(State const& s):
	S(s.S),
	A(new Matrice[S->N_spin]),
	s(new unsigned int[s.S->N_site]),
	wis(new unsigned int[s.S->N_site]),
	det(s.det)
{
	//std::cout<<"copie"<<std::endl;
	for(unsigned int i(0); i< S->N_spin; i++){
		this->A[i] = s.A[i];
	}
	for(unsigned int i(0); i< S->N_site; i++){
		this->s[i] = s.s[i];
		this->wis[i] = s.wis[i];
	}
}

State::State(unsigned int N_m, unsigned int N_spin):
	S(NULL),
	A(new Matrice[N_spin]),
	s(new unsigned int[N_spin*N_m]),
	wis(new unsigned int[N_spin*N_m]),
	det(0.0)
{
	//std::cout<<"init As ";
	init_A(N_m,N_spin);
	//std::cout<<"minimal"<<std::endl;
}

State::~State(){
	//std::cout<<"destructeur : state"<<std::endl;
	delete[] s;
	delete[] wis;
	delete[] A;
}
/*}*/

/*operator*/
/*{*/
State& State::operator=(State const& s){
	this->S = s.S;
	for(unsigned int i(0); i < S->N_site; i++){ //attention, S->N_site doit != 0
		this->s[i] = s.s[i];
		this->wis[i] = s.wis[i];
	}
	this->det = s.det;
	for(unsigned int i(0); i< S->N_spin; i++){
		this->A[i] = s.A[i];
	}
	//std::cout<<"affectation : state"<<std::endl;
	return (*this);
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

State State::swap(unsigned int a, unsigned int b) const{
	State new_s(*this);
	
	new_s.s[wis[b]] = s[wis[a]];
	new_s.s[wis[a]] = s[wis[b]];

	new_s.wis[b] = wis[a];
	new_s.wis[a] = wis[b];
	new_s.compute_matrices();
	new_s.compute_det();
	
	return new_s;
}
/*}*/

/*methods that modify the class*/
/*{*/
void State::init_A(unsigned int N_m, unsigned int N_spin){
	for(unsigned int i(0); i<N_spin; i++){
		A[i] = Matrice (N_m);
	}
}

void State::compute_det(){
	det = 1.0;
	for(unsigned int i(0); i<S->N_spin; i++){
		Lapack Ai(A[i],'G');
		det *= Ai.det();
	}

}

void State::compute_matrices(){
	unsigned int l(0);
	for(unsigned int i(0); i < S->N_spin; i++){
		for(unsigned int j(0); j < S->N_m; j++){
			l = s[i*S->N_m+j];
			for(unsigned int k(0); k < S->N_m; k++){
				A[i](k,j) = S->U(l,k);
			}
		}
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
	for(unsigned int i(0); i<S->N_spin; i++){
		for(unsigned int j(0); j<S->N_m; j++){
			std::cout<<wis[i*S->N_m+j]<<" ";
		}
	}
	std::cout<<std::endl;
	for(unsigned int i(0);i<S->N_spin;i++){
		std::cout<<std::endl;
		A[i].print();
	}
	std::cout<<det<<std::endl<<std::endl;
}
/*}*/
