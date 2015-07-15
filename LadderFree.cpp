#include "LadderFree.hpp"

LadderFree::LadderFree(
		Vector<unsigned int> const& ref, 
		unsigned int const& N, 
		unsigned int const& m, 
		unsigned int const& n, 
		Vector<unsigned int> const& M,  
		int const& bc, 
		Vector<double> const& t):
	System(ref,N,m,n,M,bc),
	Ladder<double>(8,"ladder-free-complex"),
	t_(t)
{
	if(status_==2){
		init_fermionic();
		system_info_.text("LadderFree : all colors experience the same Hamiltonian");
		filename_ += "-t";
		for(unsigned int j(0);j<t_.size();j++){
			filename_ += ((t_(j)>0)?"+":"")+my::tostring(t_(j));
		}
	}
}

/*{method needed for running*/
void LadderFree::compute_H(){
	H_.set(n_,n_,0);
	Matrix<int> nb;
	switch(t_.size()){
		case 2:
			{
				for(unsigned int i(0);i<n_;i++){
					nb = get_neighbourg(i);
					switch(i%spuc_){ 
						case 0:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(0);
								H_(i,nb(1,0)) = nb(1,1);
							}break;
						case 1:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(0);
							}break;
						case 2:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(1);
								H_(i,nb(1,0)) = nb(1,1);
							}break;
						case 3:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(0);
							}break; 
						case 4:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(0);
								H_(i,nb(1,0)) = nb(1,1);
							}break;
						case 5:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(0);
							}break;
						case 6:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(0);
								H_(i,nb(1,0)) = nb(1,1);
							}break;
						case 7:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(1);
							}break;
						default:{ std::cerr<<"void LadderFree::compute_H(unsigned int const& c) : undefined site in unit cell"<<std::endl; }break;
					}
				}
			}break;
		case 3:
			{
				for(unsigned int i(0);i<n_;i++){
					nb = get_neighbourg(i);
					switch(i%spuc_){ 
						case 0:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(0);
								H_(i,nb(1,0)) = nb(1,1);
							}break;
						case 1:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(0);
							}break;
						case 2:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(1);
								H_(i,nb(1,0)) = nb(1,1);
							}break;
						case 3:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(2);
							}break; 
						case 4:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(0);
								H_(i,nb(1,0)) = nb(1,1);
							}break;
						case 5:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(0);
							}break;
						case 6:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(2);
								H_(i,nb(1,0)) = nb(1,1);
							}break;
						case 7:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(1);
							}break;
						default:{ std::cerr<<"void LadderFree::compute_H(unsigned int const& c) : undefined site in unit cell"<<std::endl; }break;
					}
				}
			}break;
		case 5:
			{
				for(unsigned int i(0);i<n_;i++){
					nb = get_neighbourg(i);
					switch(i%spuc_){ 
						case 0:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(0);
								H_(i,nb(1,0)) = nb(1,1);
							}break;
						case 1:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(0);
							}break;
						case 2:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(1);
								H_(i,nb(1,0)) = nb(1,1);
							}break;
						case 3:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(2);
							}break; 
						case 4:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(0);
								H_(i,nb(1,0)) = nb(1,1);
							}break;
						case 5:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(0);
							}break;
						case 6:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(3);
								H_(i,nb(1,0)) = nb(1,1);
							}break;
						case 7:
							{
								H_(i,nb(0,0)) = nb(0,1)*t_(4);
							}break;
						default:{ std::cerr<<"void LadderFree::compute_H(unsigned int const& c) : undefined site in unit cell"<<std::endl; }break;
					}
				}
			}break;
		default:{  std::cerr<<"void LadderFree::compute_H(unsigned int const& c) : no wavefunction definded for "<<t_.size()-1<<" free parameters"<<std::endl; }
	}
	H_ += H_.transpose(); 
}

void LadderFree::create(){
	compute_H();
	diagonalize(true);
	if(status_==1){
		for(unsigned int c(0);c<N_;c++){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		}
	}
}

void LadderFree::save() const {
	GenericSystem<double>::save();
	std::string t_string("");
	for(unsigned int i(0);i<t_.size()-1;i++){
		t_string += my::tostring(t_(i))+",";
	}
	t_string += my::tostring(t_.back());
	jd_write_->write("t ("+t_string+")",t_);
}
/*}*/

/*{method needed for checking*/
void LadderFree::check(){
	compute_H();
	std::cout<<"liens :"<<std::endl;
	std::cout<<links_<<std::endl;

	for(unsigned int i(0);i<n_;i++){
		std::cout << "get_neighbourg. right - top/bot - left" << std::endl;
		std::cout<<"i="<<i<<std::endl;
		std::cout<<get_neighbourg(i)<<std::endl;// shows the links
	} 
	std::cout<<"Hamiltonien"<<std::endl;
	std::cout<<H_<<std::endl;

	lattice();
	plot_band_structure();
}

void LadderFree::lattice(){
	compute_H();
	Matrix<int> nb;
	std::string color("black");
	std::string linestyle("solid");
	std::string linewidth("1pt");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps("./","lattice");
	ps.add("\\begin{pspicture}(-9,-10)(16,10)%"+filename_);
	for(unsigned int i(0);i<n_;i++) {
		xy0(0) = i/2;
		xy0(1) = i%2;
		ps.put(xy0(0)-0.20,xy0(1)+0.15,my::tostring(i));
		nb = get_neighbourg(i);

		if(nb(0,1)<0){ color = "red"; } 
		else { color = "black"; }
		xy1(0) = nb(0,0)/2;
		xy1(1) = nb(0,0)%2;
		if(xy1(0)<xy0(0)){ 
			xy1(0) = xy0(0)+1;
			linestyle="dashed";
		} else{ linestyle="solid"; }
		linewidth = my::tostring(std::abs(H_(i,nb(0,0))))+"pt";
		/*x-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		if(i%2){
			color = "black";
			linestyle="solid"; 
			xy1(0) = nb(1,0)/2;
			xy1(1) = nb(1,0)%2;
			linewidth = my::tostring(std::abs(H_(i,nb(1,0))))+"pt";
			/*y-link*/ ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth="+linewidth+",linecolor="+color+",linestyle="+linestyle);
		}
	}

	Matrix<double> polygon(4,2);
	polygon(0,0)=-0.1;
	polygon(0,1)=-0.1;
	polygon(1,0)=n_/2-0.1;
	polygon(1,1)=-0.1;
	polygon(2,0)=n_/2-0.1;
	polygon(2,1)=1.1;
	polygon(3,0)=-0.1;
	polygon(3,1)=1.1;
	ps.polygon(polygon,"linecolor=green");

	polygon(0,0)=-0.1;
	polygon(0,1)=-0.1;
	polygon(1,0)=0.9;
	polygon(1,1)=-0.1;
	polygon(2,0)=0.9;
	polygon(2,1)=1.1;
	polygon(3,0)=-0.1;
	polygon(3,1)=1.1;
	ps.polygon(polygon,"linecolor=blue");

	ps.add("\\end{pspicture}");
	ps.save(true,true,true);
}
/*}*/
