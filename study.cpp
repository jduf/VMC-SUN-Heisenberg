/*!  @file study.cpp */

#include "Read.hpp"
#include "Directory.hpp"
#include "RSTFile.hpp"
#include "Gnuplot.hpp"
#include "Parseur.hpp"

void build_rst(RSTFile& rst, std::string search_in, std::string info_dir, std::string analysis_dir, std::string info_name, std::string next_analysis_dir="");
void search_jdbin(std::string path, std::string info_dir, std::string analysis_dir, RSTFile& rst);
void plot_corr(std::string filename, std::string directory, unsigned int N, unsigned int m, unsigned int n, std::string BC, double delta, Vector<double> const& corr);
void plot_long_range_corr(std::string filename, std::string directory, unsigned int N, unsigned int m, unsigned int n, std::string BC, double delta, Vector<double> const& long_range_corr);

int main(int argc, char* argv[]){
	Linux command;
	std::string root(command.pwd());
	/*info_dir must be root as the RSTFiles are saved in the parent directory*/
	std::string info_dir(root);
	std::string info_name("info");
	std::string search_in(root);
	std::string analysis_dir(root+"analysis/");

	unsigned int study;
	if(argc==1){ study = 0; }
	else {
		search_in = argv[1];
		if(search_in.find(".jdbin") == std::string::npos){ study = 1;}
		else {study = 2;}
	}

	switch(study){
		case 0:
			{
				RSTFile rst(info_dir,"README");
				Read r_readme("README");
				std::string h;
				r_readme>>h;
				rst.text(h);
				search_in += "sim/";
				build_rst(rst,search_in,info_dir,analysis_dir,info_name);
				std::cerr<<"modifer Read->eof et peut être tout ce qui va avec"<<std::endl;
			}break; 
		case 1:
			{
				if(search_in[search_in.size()-1] != '/'){ search_in += "/"; }
				info_dir += "info/";
				std::vector<std::string> tmp(string_split(search_in,'/'));
				if(tmp.size()<2){
					std::cerr<<"study : if the update of the whole sim/ directory is requested, then call '\\study' with no argument"<<std::endl;
				} else {
					for(unsigned int i(1);i<tmp.size()-1;i++){
						if(i+2==tmp.size()){/*to update the previous rst file*/
							Directory d;
							std::string tmp_local(root+"sim/");
							for(unsigned int j(1);j<=i;j++){ tmp_local += tmp[j] + "/"; }
							d.list_dir(tmp_local);
							d.sort();
							RSTFile rst(info_dir,tmp[i]);
							for(unsigned int j(0);j<d.size();j++){
								rst.hyperlink(d.get_path(j)+d.get_name(j),info_dir+tmp[i]+"/"+d.get_name(j)+".html");
								rst.nl();
							}
						}
						info_dir += tmp[i] + "/";
						analysis_dir += tmp[i] + "/";
					}
					if(tmp.size()==2){/*to update the previous REAME.rst file*/
						Directory d;
						d.list_dir(root+"sim/");
						d.sort();
						info_name = "README";
						Read r_readme(root+"README");
						std::string h;
						r_readme>>h;
						RSTFile rst(root,info_name);
						rst.text(h);
						for(unsigned int j(0);j<d.size();j++){
							rst.hyperlink(d.get_path(j)+d.get_name(j),info_dir+d.get_name(j)+".html");
							rst.nl();
						}
					}
					info_name = tmp[tmp.size()-1];
					analysis_dir += tmp[tmp.size()-1] + "/";
					search_in = root + search_in;

					RSTFile rst_next(info_dir,info_name);
					build_rst(rst_next,search_in,info_dir,analysis_dir,info_name);
				}
			}break;
		case 2:
			{
				const std::string ext(".jdbin");
				std::string filename(search_in.substr(0, search_in.size() - ext.size()));
				std::string path(root);
				info_dir += "info/";
				std::vector<std::string> tmp(string_split(filename,'/'));
				for(unsigned int i(0);i<tmp.size()-1;i++){
					path += tmp[i] + "/";
				}
				for(unsigned int i(1);i<tmp.size()-1;i++){
					info_dir += tmp[i] + "/";
				}
				filename = tmp[tmp.size()-1];

				unsigned int type;
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

				Read r(path + filename + ".jdbin");
				RSTFile rst_sim(info_dir,filename);

				r>>type>>nruns>>ref>>N>>m>>n>>bc>>param;
				for(unsigned int j(0);j<nruns;j++){
					r>>E>>DeltaE>>Nsteps>>status>>corr;
					if(type == 2){r>>long_range_corr;}
				}
				r>>E>>DeltaE>>corr>>long_range_corr;
				switch(bc){
					case -1:{BC = "A";}break;
					case 0: {BC = "O";}break;
					case 1: {BC = "P";}break;
					default:{std::cerr<<"GenericSystem : Unknown boundary condition"<<std::endl;}
				}
				plot_corr(filename,info_dir,N,m,n,BC,param,corr);
				rst_sim.figure(info_dir+filename+"-corr.png","Correlations",1000);
				if(type == 2){
					plot_long_range_corr(filename,info_dir,N,m,n,BC,param,long_range_corr);
					rst_sim.figure(info_dir+filename+"-long-range-corr.png","Long distance correlations",1000);
				}
				rst_sim.text(r.get_header());
			}break;
	}
}

void build_rst(RSTFile& rst, std::string search_in, std::string info_dir, std::string analysis_dir, std::string info_name, std::string next_analysis_dir){
	Directory d;
	d.list_dir(search_in);
	if(d.size()>0){ d.sort(); }
	info_dir += info_name + "/";
	analysis_dir += next_analysis_dir;
	Linux command;
	command("mkdir -p " + info_dir);
	command("mkdir -p " + analysis_dir);
	for(unsigned int i(0);i<d.size();i++){
		rst.hyperlink(search_in+d.get_name(i), info_dir + d.get_name(i)+".html");
		rst.nl();
		RSTFile rst_next(info_dir,d.get_name(i));
		build_rst(rst_next,d[i], info_dir, analysis_dir, d.get_name(i), d.get_name(i) + "/");
	}
	search_jdbin(search_in, info_dir, analysis_dir, rst);
}

void search_jdbin(std::string path, std::string info_dir, std::string analysis_dir, RSTFile& rst){
	Directory d;
	d.search_file_ext(".jdbin",path,false,false);
	if(d.size()>0){
		d.sort();
		std::cout<<"rep : "<<path<<std::endl;

		Write data_mean(analysis_dir+"mean.dat");
		data_mean<<"%N m n bc param E DeltaE status"<<Write::endl;
		Write data_runs(analysis_dir+"runs.dat");
		data_runs<<"%N m n bc param E DeltaE status"<<Write::endl;

		unsigned int type;
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

		Vector<double> all_param(d.size());
		Vector<std::string> all_links(d.size());

		for(unsigned int i(0); i<d.size();i++){
			std::cout<<"----->"<<d[i]<<std::endl;

			Read r(d[i]);
			RSTFile rst_sim(info_dir,d.get_name(i));

			r>>type>>nruns>>ref>>N>>m>>n>>bc>>param;
			for(unsigned int j(0);j<nruns;j++){
				r>>E>>DeltaE>>Nsteps>>status>>corr;
				if(type == 2){r>>long_range_corr;}
				data_runs<<N<<" "<<m<<" "<<n<<" "<<bc<<" "<<param<<" "<<E<<" "<<DeltaE<<" "<<status<<Write::endl;
			}
			r>>E>>DeltaE>>corr>>long_range_corr;
			data_mean<<N<<" "<<m<<" "<<n<<" "<<bc<<" "<<param<<" "<<E<<" "<<DeltaE<<" "<<status<<Write::endl;
			switch(bc){
				case -1:{BC = "A";}break;
				case 0: {BC = "O";}break;
				case 1: {BC = "P";}break;
				default:{std::cerr<<"GenericSystem : Unknown boundary condition"<<std::endl;}
			}
			plot_corr(d.get_name(i),info_dir,N,m,n,BC,param,corr);
			rst_sim.figure(info_dir+d.get_name(i)+"-corr.png","Correlations",1000);
			if(type == 2){
				plot_long_range_corr(d.get_name(i),info_dir,N,m,n,BC,param,long_range_corr);
				rst_sim.figure(info_dir+d.get_name(i)+"-long-range-corr.png","Long distance correlations",1000);
			}
			rst_sim.text(r.get_header());
			all_param(i) = param;
			all_links(i) = info_dir+d.get_name(i)+".html";
		}

		Gnuplot gp(analysis_dir,"Ep","plot");
		gp.preplot("set xlabel '$\\delta$' offset 0,1");
		gp.preplot("set ylabel '$\\dfrac{E}{n}$' rotate by 0 offset 1");
		gp.preplot("stats 'runs.dat' u 5:6 name 'row' nooutput");
		gp.preplot("stats 'runs.dat' u ($8 > 2?$5:1/0):($8 > 2?$6:1/0) name 'select' nooutput");
		gp.preplot("stats 'mean.dat' u 5:6 name 'mean' nooutput");
		gp.preplot("set xrange [row_min_x:row_max_x]");
		gp.add("'runs.dat' u 5:(($6<select_max_y && $6>select_min_y)?$6:1/0):7 w e t '$N=" +tostring(N)+ "$ $m=" +tostring(m)+ "$ $n=" +tostring(n)+ "$ bc=" + BC + "',\\\n");
		gp.add("     'mean.dat' u 5:(($6<select_max_y && $6>select_min_y)?$6:1/0):7 w e t 'averaged over "+ tostring(nruns) + "',\\\n");
		gp.add("     mean_min_y w l t sprintf('min : (%3.4f,%3.4f)',mean_pos_min_y,mean_min_y)");
		gp.save_file();
		gp.create_image(true);

		Vector<unsigned int> index(all_param.sort());
		all_links = all_links.sort(index);
		for(unsigned int i(0);i<d.size();i++){
			rst.hyperlink(tostring(all_param(i)),all_links(i));
		}
		rst.figure(analysis_dir+"Ep.png",analysis_dir+"Ep.png",1000);
	}
}

void plot_corr(std::string filename, std::string directory, unsigned int N, unsigned int m, unsigned int n, std::string BC, double delta, Vector<double> const& corr){
	Gnuplot gp(directory,filename+"-corr","plot");
	gp.preplot("set xrange [0:"+tostring(n)+"]");
	gp.preplot("set xlabel 'site' offset 0,0.5");
	gp.preplot("set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$' offset 1");
	gp.preplot("set title '$N="+tostring(N)+"$ $m="+tostring(m)+"$ $n="+tostring(n)+"$ bc="+BC+" $\\delta="+tostring(delta)+"$'");
	Vector<double> links(corr.size());
	for(unsigned int i(0);i<links.size();i++){ links(i) = 0.5+i; }
	gp.save_data(filename+"-corr.dat",links,corr);
	gp.add("notitle");
	gp.save_file();
	gp.create_image(true);
}

void plot_long_range_corr(std::string filename, std::string directory, unsigned int N, unsigned int m, unsigned int n, std::string BC, double delta, Vector<double> const& long_range_corr){
	Gnuplot gp(directory,filename+"-long-range-corr","plot");
	gp.preplot("set xlabel '$\\|i-j\\|$' offset 0,0.5");
	gp.preplot("set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(j)>$' offset 1");
	gp.preplot("set key right bottom");
	gp.preplot("set title '$N="+tostring(N)+"$ $m="+tostring(m)+"$ $n="+tostring(n)+"$ bc="+BC+" $\\delta="+tostring(delta)+"$'");
	gp.preplot("a=1.0");
	gp.preplot("b=1.0");
	gp.preplot("eta=1.0");
	gp.preplot("f(x) = "+tostring(m*m)+".0/"+tostring(N)+".0 - a/(x*x) + b*cos(2.0*pi*x*"+tostring(m)+".0/"+tostring(N)+".0)/(x**eta)");
	gp.preplot("set fit quiet");
	gp.preplot("fit f(x) '"+filename+"-long-range-corr.dat' via a,b,eta");

	gp.preplot("set xrange [0:"+tostring(long_range_corr.size()+1)+"]");
	gp.preplot("set yrange ["+tostring(1.1*long_range_corr.min())+":"+tostring(1.1*long_range_corr.max())+"]");
	Vector<double> r(long_range_corr.size());
	for(unsigned int i(0);i<r.size();i++){ r(i) = i+1; }
	gp.save_data(filename+"-long-range-corr.dat",r,long_range_corr);
	gp.add("notitle, f(x) t sprintf('$a=%3.4f, b=%3.4f, \\eta=%3.4f$',a,b,eta)");
	gp.save_file();
	gp.create_image(true);
}

