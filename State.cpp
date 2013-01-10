#include "State.hpp"

State::State(unsigned int N_spin, unsigned int N_m, arma::Mat<double> *U):
	N_spin(N_spin), 
	N_m(N_m),
	s(new unsigned int[N_m*N_spin]),
	wis(new unsigned int[N_m*N_spin]),
	det(0.0),
	A(N_spin),
	U(U)
{
	srand(time(NULL));
	unsigned int site(0);
	unsigned int N_as(N_m*N_spin);
	unsigned int available_sites[N_as];
	for(unsigned int i(0); i < N_as; i++){
		available_sites[i]  = i;	
	}

	for(unsigned int i(0); i<N_spin; i++){
		for(unsigned int j(0); j<N_m; j++){
			site = rand() % N_as;
			s[i*N_m+j] = available_sites[site];
			wis[available_sites[site]] = i*N_m+j; 	
			for(unsigned int k(site); k < N_as-1; k++){
				available_sites[k] = available_sites[k+1];
			}
			N_as--;
		}
	}
	compute_matrices();
	compute_det();
	//std::cout<<"coustructeur base"<<std::endl;
}

State::State(State const& s):
	N_spin(s.N_spin), 
	N_m(s.N_m),
	s(new unsigned int[N_m*N_spin]),
	wis(new unsigned int[N_m*N_spin]),
	det(s.det),
	A(s.A),
	U(s.U)
{
	for(unsigned int i(0); i<N_m*N_spin; i++){
		this->s[i] = s.s[i];
		this->wis[i] = s.wis[i];
	}
	//std::cout<<"coustructeur copie"<<std::endl;
}

State::State():
	N_spin(0), 
	N_m(0),
	s(NULL),
	wis(NULL),
	det(0.0),
	A(),
	U()
{
	//std::cout<<"coustructeur vide"<<std::endl;
}

State::~State(){
	delete s;
	delete wis;
	//std::cout<<"destructeur"<<std::endl;
}

State& State::operator=(State const& s){
	this->N_spin = s.N_spin; 
	this->N_m = s.N_m;
	this->s = new unsigned int[N_m*N_spin];
	this->wis = new unsigned int[N_m*N_spin];
	this->det = s.det,
	this->A = s.A;
	this->U = s.U;

	for(unsigned int i(0); i<N_m*N_spin; i++){
		this->s[i] = s.s[i];
		this->wis[i] = s.wis[i];
	}

	//std::cout<<"affectation"<<std::endl;
	return (*this);
}

State State::swap() const {
	unsigned int s1,s2,p1,p2,a,b;
	p1 = rand() % N_m;
	p2 = rand() % N_m;
	s1 = rand() % N_spin;
	s2 = rand() % N_spin;
	while(s1==s2){
		s2 = rand() % N_spin;
	}
    a = s[s1*N_m + p1];
	b = s[s2*N_m + p2];

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
	for(unsigned int i(0); i<N_spin; i++){
		det *= arma::det(A[i]);
	}
}

void State::compute_matrices(){
	for(unsigned int i(0); i<N_spin; i++){
		A[i] = arma::zeros(N_m,N_m);
		for(unsigned int j(0); j<N_m; j++){
			for(unsigned int k(0);k<N_m;k++){
				A[i](k,j) = (*U)(s[i*N_m+j],k);
			}
		}
	}
}

void State::print() const{
	std::cout<<"{";
	for(unsigned int i(0); i<N_spin; i++){
		std::cout<<"{ ";
		for(unsigned int j(0); j<N_m; j++){
			std::cout<<s[i*N_m+j]<<" ";
		}
		std::cout<<"}";
	}
	std::cout<<"}"<<std::endl;
	for(unsigned int i(0); i<N_spin; i++){
		for(unsigned int j(0); j<N_m; j++){
			std::cout<<wis[i*N_m+j]<<" ";
		}
	}
	std::cout<<std::endl;
	for(unsigned int i(0);i<N_spin;i++){
		std::cout<<std::endl;
		A[i].print();
	}
	std::cout<<det<<std::endl<<std::endl;
}

