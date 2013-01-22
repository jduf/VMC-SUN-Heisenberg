#include "State.hpp"

State::State(System *S):
	S(S),
	s(new unsigned int[S->N_site]),
	wis(new unsigned int[S->N_site]),
	det(0.0),
	A(S->N_spin)
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
	compute_matrices();
	compute_det();
	//std::cout<<"creation"<<std::endl;
}

State::State(State const& s):
	S(s.S),
	s(new unsigned int[s.S->N_site]),
	wis(new unsigned int[s.S->N_site]),
	det(s.det),
	A(s.A)
{
	for(unsigned int i(0); i< S-> N_site; i++){
		this->s[i] = s.s[i];
		this->wis[i] = s.wis[i];
	}
	//std::cout<<"copie"<<std::endl;
}

State::State(unsigned int N_site):
	S(NULL),
	s(new unsigned int[N_site]),
	wis(new unsigned int[N_site]),
	det(0.0),
	A()
{
	//std::cout<<"minimal"<<std::endl;
}

State::~State(){
	delete s;
	delete wis;
	//std::cout<<"destructeur"<<std::endl;
}

State& State::operator=(State const& s){
	this->S = s.S;
	for(unsigned int i(0); i < S->N_site; i++){ //attention, S->N_site doit != 0
		this->s[i] = s.s[i];
		this->wis[i] = s.wis[i];
	}
	this->det = s.det;
	this->A = s.A;
	//std::cout<<"affectation"<<std::endl;
	return (*this);
}

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

void State::compute_det(){
	det = 1.0;
	for(unsigned int i(0); i<S->N_spin; i++){
		det *= arma::det(A[i]);
	}
}

void State::compute_matrices(){
	for(unsigned int i(0); i < S->N_spin; i++){
		A[i] = arma::zeros(S->N_m, S->N_m);
		for(unsigned int j(0); j < S->N_m; j++){
			for(unsigned int k(0); k < S->N_m; k++){
				A[i](k,j) = S->U(s[i*S->N_m+j],k);
			}
		}
	}
}

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

