#include "YoungTableau.hpp"

YoungTableau::YoungTableau(unsigned int N):
	N(N)
{
	std::vector<unsigned int> a;
	bool add_row_(true);
	unsigned int i(0);
	unsigned int tmp(0);

	do{
		std::cout<<"number of box in the "<<++i<<" row_ : ";
		std::cin>>tmp;
		if(std::cin && tmp){
			a.push_back(tmp);
		} else {
			std::cin.clear();
			std::cin.ignore(100,'\n');
			if ( a.size() ){
				add_row_ = false;
			} else {
				std::cout<<"the tableau must contain at least one row_"<<std::endl;
				i--;
			}
		}
	} while (add_row_ && i<N);

	(*this) = YoungTableau(a,N);
}

YoungTableau::YoungTableau(std::vector<unsigned int> a, unsigned int N):
	row_(a),
	col_(a[0],0),
	N(N)
{
	for(unsigned int i(0); i<row_.size();i++){
		yt.push_back(std::vector<unsigned int> (row_[i],0));
		for(unsigned int j(0);j<row_[i];j++){
			col_[j]++;
		}
	}
}

YoungTableau::~YoungTableau(){}

void YoungTableau::reset(){
	for(unsigned int i(0);i<col_[0];i++){
		for(unsigned int j(0);j<row_[i];j++){
			yt[i][j]=0;
		}
	}
}

bool YoungTableau::is_tableau_valid(){
	if(col_[0]!=1){
		if(col_[0]>N){ return false;}
		for(unsigned int i(1);i<col_[0];i++){
			if(row_[i]>row_[i-1]){return false;}
			if(yt[i][row_[i]-1] != 0 && 
					yt[i][row_[i]-1] == yt[i-1][row_[i]-1]) {return false;}
		}
	}
	return true;
}

unsigned int YoungTableau::dimension() const{
	double out(1.0);
	for(unsigned int i(0);i<col_[0];i++){
		for(unsigned int j(0);j<row_[i];j++){
			out *= (N-i+j) / ((row_[i]-j-1.0)+(col_[j]-i-1.0)+1.0);
		}
	}
	return out;
}

std::vector<unsigned int> YoungTableau::multiplet() const{
	std::vector<unsigned int> out(N-1,0);
	if(N < col_[0]) {std::cerr<<"YoungTableau : there's an error somewhere"<<std::endl;}
	for(unsigned int i(0);i<col_[0]-1;i++){
		out[i] = row_[i] - row_[i+1];
	}
	out[col_[0]-1] = row_[col_[0]-1];
	return out;
}

std::vector<YoungTableau> YoungTableau::multiply(YoungTableau const& b, unsigned int r, unsigned int c) const{
	std::vector<YoungTableau> out;
	for(unsigned int k(r); k<col_[0]+1 && k<N; k++){
		YoungTableau tmp(add_box(k,(1+r)));
		if(tmp.is_tableau_valid()){ out.push_back(tmp); }
	}
	c++;
	if(c == b.row_[r]) {
		c=0;
		r++;
	}
	if(r<b.col_[0]) {
		std::vector<YoungTableau> out_tmp;
		for(unsigned int i(0);i<out.size();i++){
			std::vector<YoungTableau> tmp(out[i].multiply(b,r,c));
			for(unsigned int j(0); j<tmp.size(); j++){
				for(unsigned int k(0);k<out_tmp.size();k++){
					if(tmp[j] == out_tmp[k]){
						tmp.erase(tmp.begin()+j); 
					}
				}
			}
			for(unsigned int j(0); j<tmp.size(); j++){
				out_tmp.push_back(tmp[j]);
			}
		}
		return out_tmp;
	} else {
		return out;
	}
}

YoungTableau YoungTableau::add_box(unsigned int row, unsigned int b_index) const {
	std::vector<unsigned int> a(row_);
	if(row < a.size()){ a[row]++;}
	else { a.push_back(1);}
	YoungTableau tmp(a,N);
	/*! may improve this part by copying only the relevant part*/
	for(unsigned int i(0); i<col_[0] ;i++){
		for(unsigned int j(0); j<row_[i] ;j++){
			tmp.yt[i][j] = yt[i][j];
		}
	}
	tmp.yt[row][tmp.row_[row]-1]=b_index;	
	return tmp;
}

void YoungTableau::print(std::ostream& flux) const {
	for(unsigned int i(0);i<col_[0];i++){
		for(unsigned int j(0);j<row_[i];j++){
			flux<<"\u2B1C";
		}
		if(i==0){
			std::cout<<" "<<dimension()<<" (";
			std::vector<unsigned int> mult(multiplet());
			for(unsigned int i(0);i<mult.size()-1;i++){
				std::cout<<mult[i]<<",";
			}
			std::cout<<mult[mult.size()-1]<<")";
		}
		flux<<std::endl;
	}
}

bool YoungTableau::operator==(YoungTableau const& b){
	if(row_ != b.row_){ return false; }
	if(col_ != b.col_){ return false; }
	for(unsigned int i(0);i<col_[0];i++){
		for(unsigned int j(0);j<row_[i];j++){
			if(yt[i][j] != b.yt[i][j]) {return false;}
		}
	}
	return true;
}

std::vector<YoungTableau> operator*(YoungTableau const& a, YoungTableau const& b){
	std::vector<YoungTableau> out(a.multiply(b));
	for(unsigned int i(0); i<out.size(); i++){
		//out[i].test();
		out[i].reset();
	}
	return out;
}

std::vector<YoungTableau> operator*(std::vector<YoungTableau> const& list_yt, YoungTableau const& b){
	std::vector<YoungTableau> out;
	for(unsigned int i(0); i<list_yt.size();i++){
		std::vector<YoungTableau> tmp(list_yt[i].multiply(b));
		for(unsigned int j(0); j<tmp.size();j++){
			//tmp[j].test();
			out.push_back(tmp[j]);
		}
	}
	return out;
}

std::ostream& operator<<(std::ostream& flux, YoungTableau const& yt){
	yt.print(flux);
	return flux;
}

std::ostream& operator<<(std::ostream& flux, std::vector<YoungTableau> const& list_yt){
	for(unsigned int i(0); i<list_yt.size()-1; i++){
		list_yt[i].print(flux);
		flux<<" + "<<std::endl;
	}
	list_yt[list_yt.size()-1].print(flux);
	return flux;
}

void YoungTableau::test(){
	std::cout<<" |";
	for(unsigned int j(0);j<yt[0].size();j++){
		std::cout<<col_[j];
	}
	std::cout<<std::endl;
	for(unsigned int i(0);i<yt.size();i++){
		std::cout<<row_[i]<<"|";
		for(unsigned int j(0);j<yt[i].size();j++){
			if(yt[i][j]){std::cout<<yt[i][j];}
			else{std::cout<<"\u2B1C";}
		}
		if(i==0){
			std::cout<<" "<<dimension()<<" (";
			std::vector<unsigned int> mult(multiplet());
			for(unsigned int i(0);i<mult.size()-1;i++){
				std::cout<<mult[i]<<",";
			}
			std::cout<<mult[mult.size()-1]<<")";
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

//{Uncomplete code : final_check 
/*! 
  void YoungTableau::final_check(){
  std::cerr<<"YoungTableau : final_check has not been checked"<<std::endl;
  if(col_[0]>=N){ 
  unsigned int col__to_del(row_[N-1]);
  std::vector<std::vector<unsigned int> > tmp;
  for(unsigned int i(0); i<N-1; i++){
  unsigned int n_col_(row_[i]-col__to_del);
  if(n_col_){
  tmp.push_back(std::vector<unsigned int>(n_col_,0));
  for(unsigned int j(0); j<col_[0]-1;j++){
  tmp[i][j] = yt[i][j+col__to_del];
  }
  }
  }
  yt=tmp;
  } 
  }
*/
//}
