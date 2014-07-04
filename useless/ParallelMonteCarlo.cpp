
ParallelMonteCarlo::ParallelMonteCarlo(CreateSystem* CS, unsigned int nruns, unsigned int tmax):
	CS_(CS),
	tmax_(tmax),
	nruns_(nruns)
{
	E_.set_conv(true);
	//corr_.set(GS->get_n(),true);
	//if(type == 2){
		//long_range_corr_.set(GS->get_n()/3,true);
	//}
}

void ParallelMonteCarlo::run(IOFiles& w){
#pragma omp parallel for 
	for(unsigned int i=0;i<nruns_;i++){
		MonteCarlo sim(S_,tmax_);
		sim.run();
#pragma omp critical
		{
			E_.add_sample((sim.get_system())->get_energy());
			corr_.add_sample((sim.get_system())->get_corr());
			long_range_corr_.add_sample((sim.get_system())->get_long_range_corr());
			(sim.get_system())->save(w);
		}
	}

	E_.complete_analysis();
	corr_.complete_analysis();
	long_range_corr_.complete_analysis();
}

void ParallelMonteCarlo::save(IOFiles& w){
	RST rst_mean_results;
	rst_mean_results.title("Mean results (status>2)","-");
	w.add_to_header(rst_mean_results.get());
	w("E (energy per site)",E_);
	w("corr (correlation on links)",corr_);
	w("lon_range_corr (long range correlation)",long_range_corr_);
}

