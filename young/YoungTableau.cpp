#include "YoungTableau.hpp"

YoungTableau::YoungTableau(unsigned int N):
	N(N)
{
	std::vector<unsigned int> a;
	bool add_row(true);
	unsigned int i(0);
	unsigned int tmp(0);

	std::cout<<"to end the creation of a tableau, enter 0 or a non-numerical character"<<std::endl;
	do{
		std::cout<<"number of box in the "<<++i<<" row : ";
		std::cin>>tmp;
		if(std::cin && tmp){
			a.push_back(tmp);
		} else {
			std::cin.clear();
			std::cin.ignore(100,'\n');
			if ( a.size() ){
				add_row = false;
			} else {
				std::cout<<"the tableau must contain at least one row"<<std::endl;
				i--;
			}
		}
	} while (add_row);

	(*this) = YoungTableau(a,N);
}

YoungTableau::YoungTableau(std::vector<unsigned int> a, unsigned int N):
	row(a),
	col(a[0],0),
	N(N)
{
	for(unsigned int i(0); i<this->row.size();i++){
		this->yt.push_back(std::vector<unsigned int> (this->row[i],0));
		for(unsigned int j(0);j<this->row[i];j++){
			this->col[j]++;
		}
	}
}

YoungTableau::~YoungTableau(){}

void YoungTableau::reset(){
	for(unsigned int i(0);i<this->col[0];i++){
		for(unsigned int j(0);j<this->row[i];j++){
			this->yt[i][j]=0;
		}
	}
}

std::vector<YoungTableau> YoungTableau::multiply(YoungTableau const& b, unsigned int r, unsigned int c) const{
	std::vector<YoungTableau> out;
	for(unsigned int k(r); k<this->col[0]+1; k++){
		YoungTableau tmp(this->add_box(k,(1+r)));
		if(tmp.is_tableau_valid()){ out.push_back(tmp); }
	}
	c++;
	if(c == b.row[r]) {
		c=0;
		r++;
	}
	if(r<b.col[0]) {
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
	std::vector<unsigned int> a(this->row);
	if(row < a.size()){ a[row]++;}
	else { a.push_back(1);}
	YoungTableau tmp(a,this->N);
	/*! may improve this part by copying only the relevant part*/
	for(unsigned int i(0); i<this->col[0] ;i++){
		for(unsigned int j(0); j<this->row[i] ;j++){
			tmp.yt[i][j] = this->yt[i][j];
		}
	}
	tmp.yt[row][tmp.row[row]-1]=b_index;	
	return tmp;
}

bool YoungTableau::is_tableau_valid(){
	if(this->col[0]!=1){
		if(col[0]>N){return false;}
		for(unsigned int i(1);i<col[0];i++){
			if(row[i]>row[i-1]){return false;}
			if(yt[i][row[i]-1] != 0 && 
					yt[i][row[i]-1] == yt[i-1][row[i]-1]) {return false;}
		}
	}
	return true;
}

bool YoungTableau::operator==(YoungTableau const& b){
	if(this->row != b.row){ return false; }
	if(this->col != b.col){ return false; }
	for(unsigned int i(0);i<this->col[0];i++){
		for(unsigned int j(0);j<this->row[i];j++){
			if(this->yt[i][j] != b.yt[i][j]) {return false;}
		}
	}
	return true;
}

void YoungTableau::print(std::ostream& flux) const {
	for(unsigned int i(0);i<this->col[0];i++){
		for(unsigned int j(0);j<this->row[i];j++){
			flux<<"\u2B1C";
		}
		if(i==0){flux<<" "<<this->dimension();}
		flux<<std::endl;
	}
}

double YoungTableau::dimension() const{
	double out(1.0);
	for(unsigned int i(0);i<col[0];i++){
		for(unsigned int j(0);j<row[i];j++){
			out *= (N-i+j) / ((row[i]-j-1.0)+(col[j]-i-1.0)+1.0);
		}
	}
	return out;
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

//{Uncomplete code : final_check && test
/*! 
  void YoungTableau::final_check(){
  std::cerr<<"YoungTableau : final_check has not been checked"<<std::endl;
  if(this->col[0]>=N){ 
  unsigned int col_to_del(this->row[N-1]);
  std::vector<std::vector<unsigned int> > tmp;
  for(unsigned int i(0); i<N-1; i++){
  unsigned int n_col(this->row[i]-col_to_del);
  if(n_col){
  tmp.push_back(std::vector<unsigned int>(n_col,0));
  for(unsigned int j(0); j<this->col[0]-1;j++){
  tmp[i][j] = yt[i][j+col_to_del];
  }
  }
  }
  yt=tmp;
  } 
  }
  */
void YoungTableau::test(){
	std::cout<<" |";
	for(unsigned int j(0);j<yt[0].size();j++){
		std::cout<<col[j];
	}
	std::cout<<std::endl;
	for(unsigned int i(0);i<yt.size();i++){
		std::cout<<row[i]<<"|";
		for(unsigned int j(0);j<yt[i].size();j++){
			if(this->yt[i][j]){std::cout<<this->yt[i][j];}
			else{std::cout<<"\u2B1C";}
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

//}
