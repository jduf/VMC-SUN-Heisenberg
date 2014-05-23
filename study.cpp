/*!  @file study.cpp */

#include "Directory.hpp"
#include "RSTFile.hpp"
#include "Parseur.hpp"
#include "Sampling.hpp"

void build_rst(RSTFile& rst, std::string search_in, std::string info_dir, std::string analysis_dir, std::string info_name, std::string next_analysis_dir="");
double extract_jdbin(std::string info_dir, std::string path, std::string filename, IOFiles* E_file=NULL);
void search_jdbin(std::string path, std::string info_dir, std::string analysis_dir, RSTFile& rst);
void plot_corr(std::string path, std::string filename, unsigned int nruns, unsigned int length, unsigned int N, unsigned int m, unsigned int n, std::string BC, double delta);
void plot_long_range_corr(std::string path, std::string filename, unsigned int nruns, unsigned int N, unsigned int m, unsigned int n, std::string BC, double delta);

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
		if(search_in != "README"){
			if(search_in.find(".jdbin") == std::string::npos){ study = 1; }
			else { study = 2;} 
		} else {study = 3;}
	}

	switch(study){
		case 0: /*treat everything*/
			{
				RSTFile rst(info_dir,"README");
				IOFiles r_readme("README",false);
				std::string h;
				r_readme>>h;
				rst.text(h);
				search_in += "sim/";
				build_rst(rst,search_in,info_dir,analysis_dir,info_name);
				std::cerr<<"modifer Read->eof et peut être tout ce qui va avec"<<std::endl;
			}break; 
		case 1: /*treat the repository given as argument*/
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
						IOFiles r_readme(root+"README",false);
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
		case 2:  /*treat only one jdbin file*/
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

				extract_jdbin(info_dir,path,filename);
			}break;
		case 3: /*update only the README file*/
			{
				RSTFile rst(root,"README");
				IOFiles r("README",false);
				std::string h;
				r>>h;
				rst.text(h);
				Directory d;
				d.list_dir(root+"sim/");
				d.sort();
				info_dir += "info/";
				for(unsigned int j(0);j<d.size();j++){
					rst.hyperlink(d.get_path(j)+d.get_name(j),info_dir+d.get_name(j)+".html");
					rst.nl();
				}
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

		IOFiles E_file(analysis_dir+"mean.dat",true);
		E_file<<"%N m n bc param E DeltaE conv"<<IOFiles::endl;

		Vector<double> all_param(d.size());
		Vector<std::string> all_links(d.size());
		for(unsigned int i(0); i<d.size();i++){
			std::cout<<"----->"<<d[i]<<std::endl;
			all_param(i) = extract_jdbin(info_dir,d.get_path(i),d.get_name(i),&E_file);
			all_links(i) = info_dir+d.get_name(i)+".html";
		}

		Vector<unsigned int> index(all_param.sort());
		all_links = all_links.sort(index);
		for(unsigned int i(0);i<d.size();i++){
			rst.hyperlink(tostring(all_param(i)),all_links(i));
		}

		////if(d.size()>1){
			////Gnuplot gp(analysis_dir,"Ep");
			////gp+="set xlabel '$\\delta$' offset 0,1";
			////gp+="set ylabel '$\\dfrac{E}{n}$' rotate by 0 offset 1";
			////gp+="stats 'runs.dat' u 5:6 name 'row' nooutput";
			////gp+="stats 'runs.dat' u ($8 > 2?$5:1/0):($8 > 2?$6:1/0) name 'select' nooutput";
			////gp+="stats 'mean.dat' u 5:6 name 'mean' nooutput";
			////gp.xrange("row_min_x","row_max_x");
			//////gp+="plot 'runs.dat' u 5:(($6<select_max_y && $6>select_min_y)?$6:1/0):7 w e t '$N=" +tostring(N)+ "$ $m=" +tostring(m)+ "$ $n=" +tostring(n)+ "$ bc=" + BC + "',\\";
			////gp+="plot 'runs.dat' u 5:(($6<select_max_y && $6>select_min_y)?$6:1/0):7 w e t 'fix the title issue',\\";
			//////gp+="     'mean.dat' u 5:(($6<select_max_y && $6>select_min_y)?$6:1/0):7 w e t 'averaged over "+ tostring(nruns) + "',\\";
			////gp+="     'mean.dat' u 5:(($6<select_max_y && $6>select_min_y)?$6:1/0):7 w e t 'here too',\\";
			////gp+="     mean_min_y w l t sprintf('min : (%3.4f,%3.4f)',mean_pos_min_y,mean_min_y)";
			////gp.save_file();
			////gp.create_image(true);
			////rst.link_figure(analysis_dir+"Ep.png",analysis_dir+"Ep.png",analysis_dir+"Ep.gp",1000);
		////}
		if(d.size()>1){
			Gnuplot gp(analysis_dir,"Ep");
			gp+="set xlabel '$\\delta$' offset 0,1";
			gp+="set ylabel '$\\dfrac{E}{n}$' rotate by 0 offset 1";
			gp+="plot 'E.dat' u 5:(($6<select_max_y && $6>select_min_y)?$6:1/0):7 w e t 'fix the title issue'";
			//gp+="plot 'E.dat' u 5:(($6<select_max_y && $6>select_min_y)?$6:1/0):7 w e t 'fix the title issue',\\";
			//gp+="plot for [IDX=0:"+tostring(nruns-1)+"] 'E.dat' i IDX u 1:($4==1?$2:1/0):3 w errorbars lt 1 lc 3 ps 0 notitle,\\";
			//gp+="                   'E.dat' i " + tostring(nruns-1) + " IDX u 1:($4==1?$2:1/0):3 w errorbars lt 1 lc 1 ps 0 notitle,\\";
			gp.save_file();
			gp.create_image(true);
			rst.link_figure(analysis_dir+"Ep.png",analysis_dir+"Ep.png",analysis_dir+"Ep.gp",1000);
		}
	}
}

void plot_corr(std::string path, std::string filename, unsigned int nruns, unsigned int N, unsigned int m, unsigned int n, std::string BC, double delta){
	//Vector<double> poly_e(N/m,0);
	//unsigned int i(0);
	//while(i<corr.size()){
		//for(unsigned int j(0);j<N/m;j++){
			//poly_e(j) += corr(i);
			//i++;
		//}
	//}
	//poly_e /= n*m/N;
	//poly_e.sort();
	//gp+="set label '$\\Delta = "+tostring(poly_e(N/m-1)-poly_e(N/m-2))+"$' at 1,"+tostring((poly_e(N/m-1)+poly_e(N/m-2))/2.0);
	//gp+="                   '"+filename+"-corr.dat' i " + tostring(nruns) + " u 1:2:3 w errorbars lt 1 lc 2 lw 2 notitle,\\";
	//gp+=tostring(poly_e(N/m-1))+" lc 3 notitle, "+tostring(poly_e(N/m-2)) + " lc 3 notitle";

	Gnuplot gp(path,filename+"-corr");
	gp+="set xlabel 'site' offset 0,0.5";
	gp+="set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(i+1)>$' offset 1";
	gp+="set title '$N="+tostring(N)+"$ $m="+tostring(m)+"$ $n="+tostring(n)+"$ bc="+BC+" $\\delta="+tostring(delta)+"$'";
	gp+="plot for [IDX=0:"+tostring(nruns-1)+"] '"+filename+"-corr.dat' i IDX u 1:($4==1?$2:1/0):3 w errorbars lt 1 lc 3 ps 0 notitle,\\";
	gp+="     for [IDX=0:"+tostring(nruns-1)+"] '"+filename+"-corr.dat' i IDX u 1:($4==0?$2:1/0):3 w errorbars lt 1 lc 1 ps 0 notitle,\\";
	gp+="                   '"+filename+"-corr.dat' i " + tostring(nruns) + " u 1:2:3 w errorbars lt 1 lc 2 lw 2 notitle";
	gp.save_file();
	gp.create_image(true);
}

void plot_long_range_corr(std::string path, std::string filename, unsigned int nruns, unsigned int length, unsigned int N, unsigned int m, unsigned int n, std::string BC, double delta){
	Gnuplot gp(path,filename+"-long-range-corr");
	gp+="stats '"+filename+"-long-range-corr.dat' nooutput";
	gp.xrange(0,length+1);
	gp.yrange("1.1*STATS_min_y","1.1*STATS_max_y");
	gp+="set xlabel '$\\|i-j\\|$' offset 0,0.5";
	gp+="set ylabel '$<S_{\\alpha}^{\\beta}(i)S_{\\beta}^{\\alpha}(j)>$' offset 1";
	gp+="set title '$N="+tostring(N)+"$ $m="+tostring(m)+"$ $n="+tostring(n)+"$ bc="+BC+" $\\delta="+tostring(delta)+"$'";
	gp+="set key right bottom";
	gp+="a=1.0";
	gp+="b=1.0";
	gp+="eta=1.0";
	gp+="m="+tostring(m)+".0";
	gp+="N="+tostring(N)+".0";
	gp+="f(x) = a/(x*x) + b*cos(2.0*pi*x*m/N)/(x**eta)";
	gp+="set fit quiet";
	switch(N/m){
		case 2:{ gp+="fit [3:"+tostring(length)+"] f(x) '"+filename+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; } break;
		case 3:{
				   switch((length + 1) % 3){
					   case 0:{ gp+="fit [2:"+tostring(length)+"] f(x) '"+filename+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
					   case 1:{ gp+="fit [5:"+tostring(length)+"] f(x) '"+filename+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
					   case 2:{ gp+="fit [3:"+tostring(length)+"] f(x) '"+filename+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
				   }break;
			   }break;
		default :{ gp+="fit ["+tostring(N-1)+":] f(x) '"+filename+"-long-range-corr.dat' i "+tostring(nruns)+" via a,b,eta"; }break;
	}
	gp+="plot for [IDX=0:"+tostring(nruns-1)+"] '"+filename+"-long-range-corr.dat' i IDX u 1:($4==1?$2:1/0):3 w errorbars lt 1 lc 3 ps 0 notitle,\\";
	gp+="     for [IDX=0:"+tostring(nruns-1)+"] '"+filename+"-long-range-corr.dat' i IDX u 1:($4==0?$2:1/0):3 w errorbars lt 1 lc 1 ps 0 notitle,\\";
	gp+="                   '"+filename+"-long-range-corr.dat' i "+tostring(nruns)+" u 1:2:3 w errorbars lt 1 lc 2 lw 2 notitle,\\";
	gp+="                   f(x) notitle";
	gp.save_file();
	gp.create_image(true);
}

double extract_jdbin(std::string info_dir, std::string path, std::string filename, IOFiles* E_file){
	unsigned int type;
	unsigned int nruns;
	Vector<unsigned int> ref;
	unsigned int N;
	unsigned int m;
	unsigned int n;
	int bc; std::string BC;
	double param;
	Data<double> E;
	DataSet<double> corr;
	DataSet<double> long_range_corr;

	IOFiles r(path+filename+".jdbin",false);
	RSTFile rst_jdbin(info_dir,filename);
	{
		IOFiles corr_file(info_dir+filename+"-corr.dat",true);
		IOFiles long_range_corr_file(info_dir+filename+"-long-range-corr.dat",true);//should not be delcared when type!=2

		r>>type>>nruns>>ref>>N>>m>>n>>bc>>param;
		for(unsigned int i(0);i<nruns+1;i++){//!< the +1 is the averages over all runs
			r>>E>>corr>>long_range_corr;
			for(unsigned int j(0);j<corr.size();j++){
				corr_file<<j+0.5<<" "<<corr[j]<<IOFiles::endl;
			}
			corr_file<<IOFiles::endl<<IOFiles::endl;
			for(unsigned int j(0);j<long_range_corr.size();j++){
				long_range_corr_file<<j+1<<" "<<long_range_corr[j]<<IOFiles::endl;
			}
			long_range_corr_file<<IOFiles::endl<<IOFiles::endl;
			if(E_file){
				(*E_file)<<N<<" "<<m<<" "<<n<<" "<<bc<<" "<<param<<" "<<E<<IOFiles::endl;
			}
		}

		switch(bc){
			case -1:{BC = "A";}break;
			case 0: {BC = "O";}break;
			case 1: {BC = "P";}break;
			default:{std::cerr<<"GenericSystem : Unknown boundary condition"<<std::endl;}
		}
	}

	plot_corr(info_dir,filename,nruns,N,m,n,BC,param);
	rst_jdbin.link_figure(info_dir+filename+"-corr.png","Correlations",info_dir+filename+"-corr.gp",1000);
	if(type == 2){
		plot_long_range_corr(info_dir,filename,nruns,long_range_corr.size(),N,m,n,BC,param);
		rst_jdbin.link_figure(info_dir+filename+"-long-range-corr.png","Long distance correlations",info_dir+filename+"-long-range-corr.gp",1000);
	}
	rst_jdbin.text(r.get_header());
	return param;
}
