#include "VMCInterpolation.hpp"

VMCInterpolation::VMCInterpolation(VMCMinimization const& m):
	VMCMinimization(m,"INT"),
	interp_(m_->dof_)
{
	interp_.select_basis_function(7);
}

/*{public methods*/
void VMCInterpolation::init(){
	set_time();
	interp_.set_data();
	list_min_idx_.clear();

	std::cout<<"#######################"<<std::endl;
	std::string msg("VMCInterpolation");
	std::cout<<"#"<<msg<<std::endl;
	m_->info_.title(msg,'-');

	std::cout<<"#"<<get_filename()<<std::endl;
	m_->info_.item(get_filename());

	msg="contains "+my::tostring(m_->samples_.size())+" samples";
	std::cout<<"#"<<msg<<std::endl;
	m_->info_.item(msg);

	if(m_->samples_.size()){ search_minima(); }
	else {
		std::cerr<<__PRETTY_FUNCTION__<<" : empty samples_"<<std::endl;
		m_->info_.item("error : empty samples_");
	}
}

void VMCInterpolation::run(bool const& explore_around_minima){
	unsigned int lmis(list_min_idx_.size());
	if(lmis){
		std::string msg1("explore around minima");
		std::cout<<"#"<<msg1<<std::flush;
		Time chrono;

		if(explore_around_minima){
			double E(0.0);
			double dE(0.0);
			double Ee(0.0);
			double tmp(0.0);
			unsigned int j(0);
			unsigned int dir(0);
			unsigned int maxjter(10);
			unsigned int thread;
			Vector<unsigned int> pos;

			thread = omp_get_max_threads();
			RandArray<double> rnd(thread);
			RandArray<unsigned int> choose_dir(thread);
			for(unsigned int i(0);i<thread;i++){
				rnd.set(i,0.0,1.0);
				choose_dir.set(i,0,2*m_->dof_-1);
			}
#pragma omp parallel for private(E,dE,Ee,tmp,j,dir,thread,pos)
			for(unsigned int i=0;i<lmis;i++){
				pos = list_min_idx_[i];
				thread = omp_get_thread_num();
				E = 0.0;
				j = 0;
				do{
					evaluate(m_->ps_,pos,tmp,dE,Ee);
					dir = choose_dir.get(thread);
					if( tmp < E*rnd.get(thread) ){
						E = tmp;
						if(dir%2){
							if( pos(dir/2)+1<m_->ps_[dir/2].size() ){
								pos(dir/2)++;
							} else { j=maxjter; }
						} else {
							if( pos(dir/2)>1 ){
								pos(dir/2)--;
							} else { j=maxjter; }
						}
					} else { j=maxjter; }
				} while ( j++<maxjter && 2.0*j*dE/std::abs(Ee-E) > rnd.get(thread) );
			}
		} else {
#pragma omp parallel for
			for(unsigned int i=0;i<lmis;i++){
				evaluate(m_->ps_,list_min_idx_[i]);
			}
		}

		std::string msg2(" (done in "+my::tostring(chrono.elapsed())+"s)");
		std::cout<<msg2<<std::endl;
		m_->info_.item(msg1+msg2);
	}
}

void VMCInterpolation::plot(){
	if(m_->dof_<4){
		out_ = new IOFiles(get_path()+get_filename()+".dat",true,false);
		m_->samples_.set_target();
		double min(0);
		while(m_->samples_.target_next()){
			if(m_->samples_.get().get_energy().get_x()<min){ min = m_->samples_.get().get_energy().get_x(); }
			(*out_)<<m_->samples_.get().get_param()<<" "<<m_->samples_.get().get_energy()<<IOFiles::endl;
		}
		delete out_;

		out_ = new IOFiles(get_path()+get_filename()+"-spline.dat",true,false);
		Vector<unsigned int> idx(m_->dof_,0);
		while(go_through_parameter_space(m_->ps_,idx,0,0,&VMCInterpolation::save_interp_data));
		delete out_;
		out_ = NULL;

		Gnuplot plot(get_path(),get_filename());
		plot.range("x",m_->ps_[0](0),m_->ps_[0].back());
		if(m_->dof_==1){
			plot+="plot '"+get_filename()+"-spline.dat' u 1:2   lt 7 lc 1 t 'spline',\\";
			plot+="     '"+get_filename()+".dat'        u 1:2:3 lt 7 lc 2 notitle";
		} else {
			plot.range("y",m_->ps_[1](0),m_->ps_[1].back());
			plot.label("x","x");
			plot.label("y","y");
			plot.label("z","z");
			plot+="set ticslevel 0";
		}
		if(m_->dof_==2){
			plot.range("z","-1","");
			plot+="splot '"+get_filename()+"-spline.dat' u 1:2:3 lt 7 lc 1 t 'spline',\\";
			plot+="      '"+get_filename()+".dat'        u 1:2:3 lt 7 lc 2 notitle";
		}
		if(m_->dof_==3){
			plot.range("z",m_->ps_[2](0),m_->ps_[2].back());
			plot+="splot '"+get_filename()+"-spline.dat' u 1:2:3:($4>"+my::tostring(min*0.8)+"?0:exp("+my::tostring(min)+"-$4)) w p lt 7 lc 1 ps variable t 'spline',\\";
			plot+="      '"+get_filename()+".dat'        u 1:2:3:($4>"+my::tostring(min*0.8)+"?0:exp("+my::tostring(min)+"-$4)) w p lt 7 lc 2 ps variable notitle";
		}
		plot.save_file();
	} else {
		out_ = new IOFiles(get_path()+get_filename()+".dat",true,false);
		Vector<double> param(m_->dof_);
		for(unsigned int k(0);k<m_->ps_[0].size();k++){
			for(unsigned int j(0);j<m_->dof_;j++){
				param(j) = m_->ps_[j](0)-m_->ps_[j](k);
			}
			std::shared_ptr<MCSim> sim(std::make_shared<MCSim>(param));
			(*out_)<<param.norm_squared()<<" "<<interp_(param)<<" ";
			if(m_->samples_.find_in_sorted_list(sim,MCSim::sort_for_merge)){
				(*out_)<<m_->samples_.get().get_energy().get_x()<<" "<<m_->samples_.get().get_energy().get_dx()<<IOFiles::endl;
			} else {
				(*out_)<<0<<" "<<0<<IOFiles::endl;
			}
		}
		delete out_;
		out_ = NULL;

		Gnuplot plot(get_path(),get_filename());
		plot.range("x","0","");
		plot+="plot '"+get_filename()+".dat' u 1:2:3 t 'cut'";
		plot.save_file();
	}
}
/*}*/

/*{private methods*/
void VMCInterpolation::search_minima(){
	Time chrono;
	std::string msg1("computing the weights");
	std::cout<<"#"<<msg1<<std::flush;

	m_->samples_.set_target();
	while(m_->samples_.target_next()){
		interp_.add_data(m_->samples_.get().get_param(),m_->samples_.get().get_energy().get_x());
	}

	double dx(0.0);
	for(unsigned int i(0);i<m_->dof_;i++){
		double tmp(0);
		for(unsigned int j(1);j<m_->ps_[i].size();j++){
			tmp += m_->ps_[i](j)-m_->ps_[i](j-1);
		}
		dx += tmp/(m_->ps_[i].size()-1);
	}
	dx /= m_->dof_;

	if(interp_.compute_weights(dx,pow(m_->ps_size_,1./m_->dof_))){
		std::string msg2(" (done in "+my::tostring(chrono.elapsed())+"s"+(dx>1e-14?", error : "+my::tostring(dx)+")":")"));
		std::cout<<msg2<<std::endl;
		m_->info_.item(msg1+msg2);
		chrono.set();
		msg1="search for minima";
		if(m_->ps_size_ < 1e5){
			msg1="exhaustive "+msg1;
			std::cout<<"#"<<msg1<<std::flush;
#pragma omp parallel
			{
				unsigned int min0(omp_get_thread_num()*m_->ps_[0].size()/omp_get_num_threads());
				unsigned int max0((omp_get_thread_num()+1)*m_->ps_[0].size()/omp_get_num_threads());
				Vector<unsigned int> idx(m_->dof_,0);
				idx(0) = min0;
				while(go_through_parameter_space(m_->ps_,idx,min0,max0,&VMCInterpolation::select_if_min));
			}
		} else {
			msg1="random "+msg1;
			std::cout<<"#"<<msg1<<std::flush;
			RandArray<unsigned int>* rnd(new RandArray<unsigned int>[m_->dof_]);
			unsigned int tmp(omp_get_max_threads());
			for(unsigned int i(0);i<m_->dof_;i++){
				rnd[i].set(tmp);
				for(unsigned int j(0);j<tmp;j++){
					rnd[i].set(j,0,m_->ps_[i].size()-1);
				}
			}
#pragma omp parallel for
			for(unsigned int i=0;i<m_->dof_*100;i++){
				Vector<double> param(m_->dof_);
				Vector<unsigned int> pidx(m_->dof_);
				unsigned int thread(omp_get_thread_num());
				for(unsigned int j(0);j<m_->dof_;j++){
					pidx(j) = rnd[j].get(thread);
					param(j) = m_->ps_[j](pidx(j));
				}

				double tmp(interp_(param));
				double min(tmp);
				unsigned int dir;
				unsigned int const max_step(1e3);
				for(unsigned int j(0);j<max_step;j++){
					dir = 2*m_->dof_;
					for(unsigned int k(0);k<m_->dof_;k++){
						if(pidx(k)+1<m_->ps_[k].size()){
							param(k) = m_->ps_[k](pidx(k)+1);
							tmp = interp_(param);
							if(!std::isnan(tmp) && tmp<min){
								min = tmp;
								dir = 2*k;
							}
						}
						if(pidx(k)>0){
							param(k) = m_->ps_[k](pidx(k)-1);
							tmp = interp_(param);
							if(!std::isnan(tmp) && tmp<min){
								min = tmp;
								dir = 2*k+1;
							}
						}
						param(k) = m_->ps_[k](pidx(k));
					}
					if(dir != 2*m_->dof_){
						pidx(dir/2) += (dir%2?-1:1);
						param(dir/2) = m_->ps_[dir/2](pidx(dir/2));
					} else {
						j=max_step;
#pragma omp critical(VMCInterpolation__search_minima__local)
						{
							bool add_to_list_min(true);
							for(unsigned int k(0);k<list_min_idx_.size();k++){
								if(my::are_equal(pidx,list_min_idx_[k])){
									add_to_list_min = false;
									k = list_min_idx_.size();
								}
							}
							if(add_to_list_min){ list_min_idx_.push_back(pidx); }
						}
					}
				}
			}
			delete[] rnd;
			rnd = NULL;
		}

		msg2=" (found "+my::tostring(list_min_idx_.size())+" minimas in "+my::tostring(chrono.elapsed())+"s)";
		std::cout<<msg2<<std::endl;
		m_->info_.item(msg1+msg2);
	} else {
		std::string msg2(" FAIL ");
		std::cout<<msg2<<std::endl;
		m_->info_.item(msg1+msg2);
	}
}

bool VMCInterpolation::go_through_parameter_space(Vector<double>* x, Vector<unsigned int>& idx, unsigned int const& min0, unsigned int const& max0, void (VMCInterpolation::*f)(Vector<double>*, Vector<unsigned int> const&)){
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

void VMCInterpolation::select_if_min(Vector<double>* x, Vector<unsigned int> const& idx){
	double f;
	double f_tmp;
	bool is_min(true);
	Vector<double> param(m_->dof_);
	for(unsigned int i(0); i<m_->dof_;i++){ param(i) = x[i](idx(i)); }
	f = interp_(param);
	if(!std::isnan(f)){
		for(unsigned int j(0);j<m_->dof_;j++){
			if(is_min){
				if(idx(j)+1<m_->ps_[j].size()){
					param(j) = x[j](idx(j)+1);
					f_tmp = interp_(param);
					if(std::isnan(f_tmp) || f>f_tmp){ is_min = false; }
				}
				if(is_min && idx(j)>0){
					param(j) = x[j](idx(j)-1);
					f_tmp = interp_(param);
					if(std::isnan(f_tmp) || f>f_tmp){ is_min = false; }
				}
				param(j) = x[j](idx(j));
			}
		}
	} else { is_min = false; }

	if(is_min){
#pragma omp critical(VMCInterpolation__select_if_min__local)
		list_min_idx_.push_back(idx);
	}
}

void VMCInterpolation::save_interp_data(Vector<double>* x, Vector<unsigned int> const& idx){
	Vector<double> param(m_->dof_);
	for(unsigned int i(0); i<m_->dof_;i++){ param(i) = x[i](idx(i)); }
	(*out_)<<param<<" "<<interp_(param)<<IOFiles::endl;
}

void VMCInterpolation::evaluate(Vector<double>* x, Vector<unsigned int> const& idx, double& E, double& dE, double& Ee){
	Vector<double> param(m_->dof_);
	for(unsigned int i(0); i<m_->dof_;i++){ param(i) = x[i](idx(i)); }
	std::shared_ptr<MCSim> mcsim(VMCMinimization::evaluate(param,0));
	if(mcsim.get()){
		mcsim->check_conv(1e-5);
		E = mcsim->get_energy().get_x();
		dE= mcsim->get_energy().get_dx();
		Ee= interp_(param);
	} else {
		E = 0.0;
		dE= 0.0;
		Ee= 0.0;
	}
}

void VMCInterpolation::evaluate(Vector<double>* x, Vector<unsigned int> const& idx){
	Vector<double> param(m_->dof_);
	for(unsigned int i(0); i<m_->dof_;i++){ param(i) = x[i](idx(i)); }
	std::shared_ptr<MCSim> sim(VMCMinimization::evaluate(param,0));
}
/*}*/
