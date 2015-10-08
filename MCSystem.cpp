#include "MCSystem.hpp"

/*constructors and destructor*/
/*{*/
MCSystem::MCSystem(System const& S):
	System(S),
	s_(n_,m_),
	n_rnd_(0,n_-1),
	m_rnd_(0,m_-1)
{
	/*!Put the correct M_(c) particles of color c*/
	unsigned int c_tmp(0);
	Vector<unsigned int> M_tmp(M_);
	for(unsigned int p(0); p<m_; p++){
		for(unsigned int s(0); s<n_; s++){
			s_(s,p) = c_tmp;
			M_tmp(c_tmp) -= 1;
			if(!M_tmp(c_tmp)){ c_tmp++; }
		}
	}

	/*!Shuffle the particles*/
	for(unsigned int i(0);i<N_*(n_*m_)*(n_*m_);i++){
		swap();
		s_(new_s_[0],new_p_[0]) = new_c_[1];
		s_(new_s_[1],new_p_[1]) = new_c_[0];
	}
}

MCSystem::MCSystem(MCSystem const& mcsim):
	System(mcsim),
	s_(mcsim.s_),
	n_rnd_(0,n_-1),
	m_rnd_(0,m_-1)
{}

MCSystem::MCSystem(IOFiles& r):
	System(r),
	s_(r),
	n_rnd_(0,n_-1),
	m_rnd_(0,m_-1)
{}
/*}*/

/*public method*/
/*{*/
void MCSystem::swap(){
	new_s_[0] = n_rnd_.get();
	new_p_[0] = m_rnd_.get();
	new_c_[0] = s_(new_s_[0],new_p_[0]);
	do {
		new_s_[1] = n_rnd_.get();
		new_p_[1] = m_rnd_.get();
		new_c_[1] = s_(new_s_[1],new_p_[1]);
	} while(is_new_state_forbidden() || new_c_[0] == new_c_[1]);
}

void MCSystem::swap(unsigned int const& s0, unsigned int const& s1, unsigned int const& p0, unsigned int const& p1){
	new_s_[0] = s0;
	new_p_[0] = p0;
	new_c_[0] = s_(s0,p0);
	new_s_[1] = s1;
	new_p_[1] = p1;
	new_c_[1] = s_(s1,p1);
}

void MCSystem::update(){
	s_(new_s_[0],new_p_[0]) = new_c_[1];
	s_(new_s_[1],new_p_[1]) = new_c_[0];
}

void MCSystem::measure_new_step(){
	E_.set_x(0.0);
	if(!corr_types_[0].size()){/*!compute energy*/
		for(unsigned int l(0);l<link_types_[0].row();l++){
			for(unsigned int p0(0);p0<m_;p0++){
				for(unsigned int p1(0);p1<m_;p1++){
					swap(link_types_[0](l,0),link_types_[0](l,1),p0,p1);
					/*!if the new state is forbidden, r=0 and therefore there is no
					 * need to complete the else condition*/
					if(!is_new_state_forbidden()){
						E_.add(J_(l)*ratio(false));
					}
				}
			}
		}
	} else {/*!compute energy and the correlation*/
		double r;
		for(unsigned int l(0);l<link_types_[0].row();l++){
			corr_types_[0][l].set_x(0.0);
			for(unsigned int p0(0);p0<m_;p0++){
				for(unsigned int p1(0);p1<m_;p1++){
					swap(link_types_[0](l,0),link_types_[0](l,1),p0,p1);
					/*!if the new state is forbidden, r=0 and therefore there is no
					 * need to complete the else condition*/
					if(!is_new_state_forbidden()){
						r = ratio(false);
						E_.add(J_(l)*r);
						corr_types_[0][l].add(r);
					}
				}
			}
		}
	}
	E_.divide(n_);

	/*!compute the long range correlation*/
	double diag_term(1.0*m_*m_/N_);
	for(unsigned int i(1);i<n_corr_;i++){
		for(unsigned int l(0);l<corr_types_[i].size();l++){
			corr_types_[i][l].set_x(-diag_term);
			for(unsigned int p0(0);p0<m_;p0++){
				for(unsigned int p1(0);p1<m_;p1++){
					swap(link_types_[i](l,0),link_types_[i](l,1),p0,p1);
					if(!is_new_state_forbidden() && new_c_[0] == new_c_[1]){ corr_types_[i][l].add(1.0); }
				}
			}
		}
	}
}

void MCSystem::add_sample(){
	E_.add_sample();
	for(unsigned int i(0);i<n_corr_;i++){
		corr_types_[i].add_sample();
	}
}

void MCSystem::write(IOFiles& w) const {
	w<<s_;
}
/*}*/

/*private methods*/
/*{*/
bool MCSystem::is_new_state_forbidden(){
	for(unsigned int i(0); i<m_; i++){
		if(i != new_p_[0] && s_(new_s_[0],i) == new_c_[1]){ return true; }
		if(i != new_p_[1] && s_(new_s_[1],i) == new_c_[0]){ return true; }
	}
	return false;
}
/*}*/
