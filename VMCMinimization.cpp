#include "VMCMinimization.hpp"

VMCMinimization::VMCMinimization(Parseur& P):
	tmax_(P.get<unsigned int>("tmax")),
	basename_("")
{
	unsigned int i(0);
	IOFiles* in(P.find("load",i,false)?(new IOFiles(P.get<std::string>(i),false)):NULL);
	std::string wf(in?in->read<std::string>():P.get<std::string>("wf"));
	unsigned int N(in?in->read<unsigned int>():P.get<unsigned int>("N"));
	unsigned int m(in?in->read<unsigned int>():P.get<unsigned int>("m"));
	unsigned int n(in?in->read<unsigned int>():P.get<unsigned int>("n"));
	int bc(in?in->read<int>():P.get<int>("bc"));

	system_param_.set("wf",wf);
	system_param_.set("N",N);
	system_param_.set("m",m);
	system_param_.set("n",n);
	system_param_.set("bc",bc);

	basename_ += "-wf"+system_param_.get<std::string>("wf");
	basename_ += "-N" +my::tostring(N);
	basename_ += "-m" +my::tostring(m);
	basename_ += "-n" +my::tostring(n);
	basename_ += "-bc"+my::tostring(bc);

	if(in){
		unsigned int size(in->read<int>());
		while(size--){ all_results_.add_end(std::make_shared<MCSim>(*in)); }
		pso_info_.text("loads data from "+in->get_filename()+RST::nl_);

		delete in;
		in = NULL;
	}
}

void VMCMinimization::complete_analysis(double const& converged_criterion){
	all_results_.set_target();
	while ( all_results_.target_next() ){
		all_results_.get().complete_analysis(converged_criterion);
	}
}

void VMCMinimization::save() const {
	IOFiles out(get_filename()+".jdbin",true);
	out.write("wf",system_param_.get<std::string>("wf"));
	out.write("N", system_param_.get<unsigned int>("N"));
	out.write("m", system_param_.get<unsigned int>("m"));
	out.write("n", system_param_.get<unsigned int>("n"));
	out.write("bc",system_param_.get<int>("bc"));
	out.write("#", all_results_.size());
	out.add_header()->nl();
	out.add_header()->text(pso_info_.get());
	while(all_results_.target_next()){ all_results_.get().write(out); }
}
