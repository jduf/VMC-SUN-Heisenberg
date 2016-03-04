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
			M_tmp(c_tmp)--;
			if(!M_tmp(c_tmp)){ c_tmp++; }
		}
	}

	/*!Shuffle the particles*/
	for(unsigned int i(0);i<N_*(n_*m_)*(n_*m_);i++){
		swap();
		s_(new_s_[0],new_p_[0]) = new_c_[1];
		s_(new_s_[1],new_p_[1]) = new_c_[0];
	}

	set_list_of_measures();
}

MCSystem::MCSystem(MCSystem const& mcsim):
	System(mcsim),
	s_(mcsim.s_),
	n_rnd_(0,n_-1),
	m_rnd_(0,m_-1),
	list_of_measures_(mcsim.list_of_measures_)
{}

MCSystem::MCSystem(IOFiles& r):
	System(r),
	s_(r),
	n_rnd_(0,n_-1),
	m_rnd_(0,m_-1)
{ set_list_of_measures(); }
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
	} while (is_new_state_forbidden() || new_c_[0] == new_c_[1]);
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
	unsigned int L(obs_[0].nlinks());
	int* idx[3];
	idx[0] = obs_[0].get_links().ptr();
	idx[1] = idx[0]+L;
	idx[2] = idx[1]+L;

	Observable* O;
	Data<double>* E(&obs_[0][0]);
	E->set_x(0.0);
	for(unsigned int i(0);i<list_of_measures_.size();i++){
		switch(list_of_measures_[i].first){
			case 0:// only for energy
				{
					if(m_==1){
						for(unsigned int l(0);l<L;l++){
							swap(*idx[0],*idx[1],0,0);
							E->add(J_(l)*ratio());
							idx[0]++;
							idx[1]++;
						}
					} else {
						for(unsigned int l(0);l<L;l++){
							for(unsigned int p0(0);p0<m_;p0++){
								for(unsigned int p1(0);p1<m_;p1++){
									swap(*idx[0],*idx[1],p0,p1);
									/*!if the new state is forbidden, r=0 and therefore
									 * there is no need to complete the else condition*/
									if(!is_new_state_forbidden()){ E->add(J_(l)*ratio()); }
								}
							}
							idx[0]++;
							idx[1]++;
						}
					}
					E->divide(n_);
				}break;
			case 1:// for energy and bond energy
				{
					double r;
					O = &obs_[list_of_measures_[i].second];
					O->set_x(0.0);
					if(m_==1){
						for(unsigned int l(0);l<L;l++){
							swap(*idx[0],*idx[1],0,0);
							r = J_(l)*ratio();
							E->add(r);
							O->add(*idx[2],r);
							idx[0]++;
							idx[1]++;
							idx[2]++;
						}
					} else {
						for(unsigned int l(0);l<L;l++){
							for(unsigned int p0(0);p0<m_;p0++){
								for(unsigned int p1(0);p1<m_;p1++){
									swap(*idx[0],*idx[1],p0,p1);
									/*!if the new state is forbidden, r=0 and therefore
									 * there is no need to complete the else condition*/
									if(!is_new_state_forbidden()){
										r = J_(l)*ratio();
										E->add(r);
										O->add(*idx[2],r);
									}
								}
							}
							idx[0]++;
							idx[1]++;
							idx[2]++;
						}
					}
					E->divide(n_);
				}break;
			case 2:// for long range correlation 
				{
					O = &obs_[list_of_measures_[i].second];
					O->set_x(-1.0*m_*m_/N_);
					L = O->nlinks();
					idx[0] = O->get_links().ptr();
					idx[1] = idx[0]+L;
					idx[2] = idx[1]+L;

					if(m_ == 1){
						for(unsigned int l(0);l<L;l++){
							if(s_(*idx[0],0) == s_(*idx[1],0)){ O->add(*idx[2],1.0); }
							idx[0]++;
							idx[1]++;
							idx[2]++;
						}
					} else {
						for(unsigned int l(0);l<L;l++){
							for(unsigned int p0(0);p0<m_;p0++){
								for(unsigned int p1(0);p1<m_;p1++){
									if(s_(*idx[0],p0) == s_(*idx[1],p1)){ O->add(*idx[2],1.0); }
								}
							}
							idx[0]++;
							idx[1]++;
							idx[2]++;
						}
					}
				}break;
			case 3:// for color occupation
				{
					O = &obs_[list_of_measures_[i].second]; 
					O->set_x(0); 
					if(m_ == 1){
						for(unsigned int i(0);i<n_;i++){
							O->add((*O)(i,s_(i,0)),1.0);
						}
					} else {
						for(unsigned int i(0);i<n_;i++){
							for(unsigned int p(0);p<m_;p++){
								O->add((*O)(i,s_(i,p)),1.0);
							}
						}
						std::cerr<<__PRETTY_FUNCTION__<<" : need to check"<<std::endl;
					}
				}break;
		}
	}
}

void MCSystem::add_sample(){
	for(unsigned int i(0);i<obs_.size();i++){ obs_[i].add_sample(); }
}

void MCSystem::write(IOFiles& w) const {
	w<<s_;
}
/*}*/

/*private methods*/
/*{*/
bool MCSystem::is_new_state_forbidden(){
	for(unsigned int i(0);i<m_;i++){
		if(i != new_p_[0] && s_(new_s_[0],i) == new_c_[1]){ return true; }
		if(i != new_p_[1] && s_(new_s_[1],i) == new_c_[0]){ return true; }
	}
	return false;
}

void MCSystem::set_list_of_measures(){
	list_of_measures_.push_back(std::pair<unsigned int,unsigned int>(0,0));
	for(unsigned int i(0);i<obs_.size();i++){
		switch(obs_[i].get_type()){
			case 1://set bond energy
				{ list_of_measures_[0] = std::pair<unsigned int,unsigned int>(1,i); }break;
			case 2://set long range correlations
				{ list_of_measures_.push_back(std::pair<unsigned int,unsigned int>(2,i)); }break;
			case 3://set color occupation
				{ list_of_measures_.push_back(std::pair<unsigned int,unsigned int>(3,i)); }break;
		}
	}
}
/*}*/
