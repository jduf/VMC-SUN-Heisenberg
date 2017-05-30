#include "VMCSystematic.hpp"

VMCSystematic::VMCSystematic(VMCMinimization const& m):
	VMCMinimization(m,"SYS")
{}

VMCSystematic::VMCSystematic(IOFiles& in):
	VMCMinimization(in,true,"SYS")
{}

/*{public methods*/
void VMCSystematic::run(double const& dEoE, unsigned int const& tmax, unsigned int const& maxiter, std::string const& save_in){
	tmax_ = tmax;
	dEoE_ = dEoE;
	if(tmax_ && dEoE_){
		std::cout<<"#######################"<<std::endl;
		std::string msg("VMCSystematic");
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.title(msg,'-');

		msg = "do a systematic measure over the phase space";
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);
		msg = my::tostring(m_->ps_size_)+" samples (max time "+my::tostring(tmax_*maxiter*m_->ps_size_)+"s)";
		std::cout<<"#"<<msg<<std::endl;
		m_->info_.item(msg);

		maxiter_ = maxiter;
		Vector<unsigned int> idx(m_->dof_,0);
		total_eval_ = m_->ps_size_;
		progress_ = 0;
		while(go_through_parameter_space(m_->ps_,idx,0,0,&VMCSystematic::evaluate));
		save(save_in);
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : tmax_ = 0"<<std::endl; }
}

void VMCSystematic::plot(){
	if(m_->samples_.size()){
		IOFiles data(get_path()+get_filename()+".dat",true,false);
		m_->samples_.set_target();
		Vector<double> param;
		while(m_->samples_.target_next()){
			param = m_->samples_.get().get_param();
			data<<m_->samples_.get().get_param()<<" "<<m_->samples_.get().get_energy()<<IOFiles::endl;
		}

		switch(m_->dof_){
			case 1:
				{
					Gnuplot gp(get_path(),get_filename());
					gp.label("y2","$\\frac{E}{n}$","rotate by 0");
					gp+="plot '"+get_filename()+".dat' u 1:2:3 w e notitle";
					gp.save_file();
					gp.create_image(true,"png");
				}break;
			case 2:
				{
					Gnuplot gp(get_path(),get_filename());
					gp.label("z","$\\frac{E}{n}$");
					gp+="splot '"+get_filename()+".dat' u 1:2:($3+$4) lc 1 pt 8 notitle,\\";
					gp+="      '"+get_filename()+".dat' u 1:2:($3-$4) lc 1 pt 7 notitle";
					gp.save_file();
					gp.create_image(true,"png");
				}break;
			default:
				{
					Gnuplot gp(get_path(),get_filename());
					gp+="plot '"+get_filename()+".dat' u  w e notitle";
					gp.save_file();
				}break;
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : no samples"<<std::endl; }
}

void VMCSystematic::analyse(std::string const& path, std::string const& filename, List<MCSim>& kept_samples) const {
	if(m_->samples_.size()){
		Vector<unsigned int> ref(m_->s_->get_ref());
		switch(ref(0)){
			case 2:
				{
					IOFiles data(path+filename+".dat",true,false);

					List<MCSim>::Node* min(NULL);
					m_->samples_.target_next();
					while(m_->samples_.target_next()){
						data<<m_->samples_.get().get_param()<<" "<<m_->samples_.get().get_energy()<<IOFiles::endl;

						if(!min || min->get()->get_energy().get_x()>m_->samples_.get().get_energy().get_x())
						{ min = m_->samples_.get_target(); }
					}
					if(min){ kept_samples.add_end(min->get()); }
					else { std::cerr<<__PRETTY_FUNCTION__<<" : no minium found for "<<path+filename+".jdbin"<<std::endl; }
				}break;
			case 3:
				{
					IOFiles data(path+filename+".dat",true,false);

					List<MCSim>::Node* min(NULL);
					m_->samples_.target_next();
					while(m_->samples_.target_next()){
						data<<m_->samples_.get().get_param()<<" "<<m_->samples_.get().get_energy()<<IOFiles::endl;

						if(
								m_->samples_.prev() &&
								m_->samples_.next() &&
								m_->samples_.prev()->get()->get_energy().get_x() > m_->samples_.get().get_energy().get_x() &&
								m_->samples_.next()->get()->get_energy().get_x() > m_->samples_.get().get_energy().get_x()
						  )
						{
							if(!min || min->get()->get_energy().get_x()>m_->samples_.get().get_energy().get_x())
							{ min = m_->samples_.get_target(); }
						}
					}
					if(min){ kept_samples.add_end(min->get()); }

					Gnuplot gp(path,filename);
					gp+="plot '"+filename+".dat' u 1:2:3 w e notitle";
					gp.save_file();
					gp.create_image(true,"png");
				}break;
			case 6:
				{
					if(ref(1) == 1 && (ref(2) == 3 || ref(2) == 4)){
						IOFiles data(path+filename+".dat",true,false);

						/*keep the two minima for at (0pp,000) and (p00,pp) for the honeycomb wave functions */
						List<MCSim>::Node* min_neg(NULL);
						List<MCSim>::Node* min_pos(NULL);
						m_->samples_.target_next();
						while(m_->samples_.target_next()){
							data<<m_->samples_.get().get_param()<<" "<<m_->samples_.get().get_energy()<<IOFiles::endl;

							if(
									m_->samples_.prev() &&
									m_->samples_.next() &&
									m_->samples_.prev()->get()->get_energy().get_x() > m_->samples_.get().get_energy().get_x() &&
									m_->samples_.next()->get()->get_energy().get_x() > m_->samples_.get().get_energy().get_x()
							  )
							{
								if(m_->samples_.get().get_param()(0)<0){
									if(!min_neg || min_neg->get()->get_energy().get_x()>m_->samples_.get().get_energy().get_x())
									{ min_neg = m_->samples_.get_target(); }
								}
								if(m_->samples_.get().get_param()(0)>0){
									if(!min_pos || min_pos->get()->get_energy().get_x()>m_->samples_.get().get_energy().get_x())
									{ min_pos = m_->samples_.get_target(); }
								}
							}
						}
						if(min_neg){ kept_samples.add_end(min_neg->get()); }
						if(min_pos){ kept_samples.add_end(min_pos->get()); }

						Gnuplot gp(path,filename);
						gp+="plot '"+filename+".dat' u 1:2:3 w e notitle";
						gp.save_file();
						gp.create_image(true,"png");
					} else { std::cerr<<__PRETTY_FUNCTION__<<" : undefined analyse for "<<ref<<std::endl; }
				}break;
			default:
				{ std::cerr<<__PRETTY_FUNCTION__<<" : undefined analyse for "<<ref<<std::endl; }
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : no samples"<<std::endl; }
}
/*}*/

/*{private methods*/
bool VMCSystematic::go_through_parameter_space(Vector<double>* x, Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0, void (VMCSystematic::*f)(Vector<double>*, Vector<unsigned int> const&)){
	(this->*f)(x,idx);

	idx(0)++;
	if(min0==0 && max0==0){
		for(unsigned int i(0);i<m_->dof_;i++){
			if(idx(i) == x[i].size()){
				if(i+1 == m_->dof_){ return false; }
				idx(i) = 0;
				idx(i+1)++;
			}
		}
	} else {
		if(idx(0) == max0){
			if(1 == m_->dof_){ return false; }
			idx(0) = min0;
			idx(1)++;
		}
		for(unsigned int i(1);i<m_->dof_;i++){
			if(idx(i) == x[i].size()){
				if(i+1 == m_->dof_){ return false; }
				idx(i) = 0;
				idx(i+1)++;
			}
		}
	}

	return go_through_parameter_space(x,idx,min0,max0,f);
}

void VMCSystematic::evaluate(Vector<double>* x, Vector<unsigned int> const& idx){
	Vector<double> param(m_->dof_);
	for(unsigned int i(0); i<m_->dof_;i++){ param(i) = x[i](idx(i)); }
	evaluate_until_precision(param,0,dEoE_,tmax_,maxiter_);
}
/*}*/
