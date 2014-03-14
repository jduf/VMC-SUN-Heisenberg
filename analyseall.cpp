/*!  @file analyseall.cpp */

#include "Read.hpp"
#include "Directory.hpp"
#include "CreateSystem.hpp"
#include "RSTFile.hpp"
#include "Gnuplot.hpp"

void analyse_monte_carlo(std::string sim_name, std::string sim_directory, std::string info_directory, std::string analyse_directory);
void analyse_monte_carlo_parameters(std::string sim_name, std::string sim_directory, std::string info_directory, std::string analyse_directory);
void analyse_corr(std::string directory, std::string filename, unsigned int N, unsigned int m, unsigned int n, std::string BC, double delta, Vector<double> const& corr);
void analyse_long_range_corr(std::string directory, std::string filename, unsigned int N, unsigned int m, unsigned int n, std::string BC, double delta, Vector<double> const& long_range_corr);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	std::string root(P.get<std::string>("0"));
	if(!P.status()){
		Linux command;
		if(root == "."){ root = command.pwd(); }
		else {root = command.pwd()+root; }
		if(root[root.size()-1] != '/'){ root += "/"; }
		std::string save_in("info/");
		command("mkdir -p " + root + save_in);

		RSTFile readme("README",root+save_in);
		Read r_readme("README");
		std::string h;
		r_readme>>h;
		readme.text(h);

		Directory d;
		d.search_file_ext("-mc.jdbin",root+"sim/",false,false);
		d.sort();

		std::string relative_path;
		unsigned int type;
		for(unsigned int i(0); i<d.size();i++){
			std::cout<<"Simulation : "<<d[i]<<std::endl;
			Read r(d[i]);
			r>>type>>relative_path;
			switch(type){
				case 1:{
						   analyse_monte_carlo_parameters(d.get_name(i),root+"sim/"+relative_path,root+"info/",root+"analysis/");
						   readme.hyperlink(d.get_name(i),root+"info/"+d.get_name(i)+".html");
						   readme.nl();
					   }break;
				case 2:{ 
						   std::string sim_name;
						   r>>sim_name;
						   analyse_monte_carlo(sim_name,root+"sim/"+relative_path,root+"info/",root+"analysis/");
						   readme.hyperlink(sim_name,root+"info/"+sim_name+".html");
						   readme.nl();
					   }break;
				default:{std::cout<<"unkown simulation"<<std::endl;}
			}
		}

		std::cerr<<"modifer Read->eof et peut être tout ce qui va avec"<<std::endl;
	}
}

void analyse_monte_carlo(std::string sim_name, std::string sim_directory, std::string info_directory, std::string analyse_directory){
	RSTFile rst_sim(sim_name,info_directory);
	unsigned int nruns;
	unsigned int N;
	unsigned int m;
	unsigned int n;
	int bc;
	std::string BC;
	double param;
	double E;
	double DeltaE;
	unsigned int Nsteps;
	unsigned int status;
	Vector<double> corr;
	Vector<double> long_range_corr;
	Vector<unsigned int> ref;

	Write data_runs(analyse_directory+sim_name+"-runs.dat");
	Read r(sim_directory+sim_name+".jdbin");
	r>>nruns>>ref>>N>>m>>n>>bc>>param;
	for(unsigned int j(0);j<nruns;j++){
		r>>E>>DeltaE>>Nsteps>>status>>corr>>long_range_corr;
		data_runs<<N<<" "<<m<<" "<<n<<" "<<bc<<" "<<param<<" "<<E<<" "<<DeltaE<<" "<<status<<Write::endl;
	}
	switch(bc){
		case -1:{BC = "A"; }break;
		case 0: {BC = "O"; }break;
		case 1: {BC = "P"; }break;
		default:{std::cerr<<"GenericSystem : Unknown boundary condition"<<std::endl;}
	}
	r>>E>>DeltaE>>corr>>long_range_corr;
	analyse_long_range_corr(sim_name,info_directory,N,m,n,BC,param,long_range_corr);

	std::string h(r.get_header());
	std::string data("");

	RSTFile rst_run(sim_name,info_directory);
	rst_run.text(h);

	rst_sim.figure(sim_name+".png","Long distance correlations",1000);
	rst_sim.hyperlink(tostring(param),info_directory+sim_name+".html");
}

void analyse_monte_carlo_parameters(std::string sim_name, std::string sim_directory, std::string info_directory, std::string analyse_directory){
	RSTFile rst_sim(sim_name,info_directory);
	Directory d;
	d.search_file_ext(".jdbin",sim_directory,false,false);
	d.sort();

	unsigned int nruns;
	unsigned int N;
	unsigned int m;
	unsigned int n;
	int bc;
	std::string BC;
	double param;
	double E;
	double DeltaE;
	unsigned int Nsteps;
	unsigned int status;
	Vector<double> corr;
	Vector<double> long_range_corr;
	Vector<unsigned int> ref;

	Write data_runs(analyse_directory+sim_name+"-runs.dat");
	Write data_mean(analyse_directory+sim_name+"-mean.dat");
	data_runs<<"%N m n bc param E DeltaE status"<<Write::endl;
	data_mean<<"%N m n bc param E DeltaE status"<<Write::endl;
	for(unsigned int i(0); i<d.size();i++){
		std::cout<<"-----------> "<<d[i]<<std::endl;
		Read r(d[i]);
		r>>nruns>>ref>>N>>m>>n>>bc>>param;
		for(unsigned int j(0);j<nruns;j++){
			r>>E>>DeltaE>>Nsteps>>status>>corr>>long_range_corr;
			data_runs<<N<<" "<<m<<" "<<n<<" "<<bc<<" "<<param<<" "<<E<<" "<<DeltaE<<" "<<status<<Write::endl;
		}
		switch(bc){
			case -1:{BC = "A"; }break;
			case 0: {BC = "O"; }break;
			case 1: {BC = "P"; }break;
			default:{std::cerr<<"GenericSystem : Unknown boundary condition"<<std::endl;}
		}
		r>>E>>DeltaE>>corr;
		data_mean<<N<<" "<<m<<" "<<n<<" "<<bc<<" "<<param<<" "<<E<<" "<<DeltaE<<" "<<status<<Write::endl;
		//CreateSystem CS(N,n,m,bc,param,ref);
		//CS.study(E,DeltaE,corr,info_directory);
		analyse_corr(d.get_name(i),info_directory,N,m,n,BC,param,corr);

		std::string h(r.get_header());
		std::string data("");

		RSTFile rst_run(d.get_name(i),info_directory);
		rst_run.text(h);
		rst_run.figure(d.get_name(i)+".png","Correlations",1000);
		rst_run.textit(d[i]);
		rst_sim.hyperlink(tostring(param),info_directory+d.get_name(i)+".html");
	}
	Gnuplot gp(analyse_directory,sim_name,"plot",true);
	gp.preplot("set xlabel '$\\delta$' offset 0,1");
	gp.preplot("set ylabel '$\\dfrac{E}{n}$' rotate by 0 offset 1");
	gp.preplot("stats '"+ sim_name+"-runs.dat' u 5:6 name 'row' nooutput");
	gp.preplot("stats '"+ sim_name+"-runs.dat' u ($8 > 2?$5:1/0):($8 > 2?$6:1/0) name 'select' nooutput");
	gp.preplot("set xrange [row_min_x:row_max_x]");
	gp.add("'"+sim_name+"-runs.dat' u 5:(($6<select_max_y && $6>select_min_y)?$6:1/0):7 w e t '$N=" +tostring(N)+ "$ $m=" +tostring(m)+ "$ $n=" +tostring(n)+ "$ bc=" + BC + "',\\\n");
	gp.add("     '"+sim_name+"-mean.dat' u 5:(($6<select_max_y && $6>select_min_y)?$6:1/0):7 w e t 'averaged over "+ tostring(nruns) + "',\\\n");
	gp.add("     select_min_y w l t sprintf(\"min : %3.4f\",select_min_y)");
	rst_sim.figure(analyse_directory+sim_name+".png",sim_name,1000);
}

void analyse_corr(std::string filename, std::string directory, unsigned int N, unsigned int m, unsigned int n, std::string BC, double delta, Vector<double> const& corr){
	Gnuplot gp(directory,filename,"plot",true);
	gp.preplot("set xrange [0:"+tostring(n)+"]");
	gp.preplot("set xlabel 'site' offset 0,0.5");
	gp.preplot("set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$' offset 1");
	gp.preplot("set title '$N="+tostring(N)+"$ $m="+tostring(m)+"$ $n="+tostring(n)+"$ bc="+BC+" $\\delta="+tostring(delta)+"$'");
	Vector<double> links(corr.size());
	for(unsigned int i(0);i<links.size();i++){ links(i) = 0.5+i; }
	gp.save_data(filename+"-corr.dat",links,corr);
	gp.add("notitle");
}

void analyse_long_range_corr(std::string filename, std::string directory, unsigned int N, unsigned int m, unsigned int n, std::string BC, double delta, Vector<double> const& long_range_corr){
	Gnuplot gp(directory,filename,"plot",false);
	//gp.preplot("set xrange [0:"+tostring(n)+"]");
	gp.preplot("set xlabel '$\\|i-j\\|$' offset 0,0.5");
	gp.preplot("set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(j)>$' offset 1");
	gp.preplot("set key right bottom");
	gp.preplot("set title '$N="+tostring(N)+"$ $m="+tostring(m)+"$ $n="+tostring(n)+"$ bc="+BC+" $\\delta="+tostring(delta)+"$'");
	gp.preplot("a=1.0");
	gp.preplot("b=1.0");
	gp.preplot("eta = 2.0*(1.0-1.0/"+tostring(N)+")");
	gp.preplot("f(x) = "+tostring(1.0/N)+" - a/(x*x) + b*cos(2.0*pi*x/4.0)/(x**eta)");
	gp.preplot("fit f(x) 'chain-polymerized-N4-m1-n100-P-delta0-long-range-corr.dat' via a,b");
	Vector<double> r(long_range_corr.size());
	for(unsigned int i(0);i<r.size();i++){ r(i) = i+1; }
	gp.save_data(filename+"-long-range-corr.dat",r,long_range_corr);
	gp.add("notitle, f(x) t sprintf('fit a=%3.3f, b=%3.3f',a,b)");
}
