#include "YoungTableau.hpp"

YoungTableau::YoungTableau(std::vector<unsigned int> a):
	boxes_in_row(a),
	boxes_in_col(a[0],0),
	n_box(0),
	N(3)
{
	for(unsigned int i(0); i<a.size();i++){
		yt.push_back(std::vector<unsigned int> (a[i],0));
		for(unsigned int j(0);j<a[i];j++){
			boxes_in_col[j]++;
		}
	}
	for(unsigned int i(0);i<boxes_in_row.size();i++){
		std::cout<<boxes_in_row[i]<<std::endl;
	}
	std::cout<<std::endl;
	for(unsigned int i(0);i<boxes_in_col.size();i++){
		std::cout<<boxes_in_col[i]<<std::endl;
	}
	std::cout<<std::endl;
	if(this->valid()){
		std::cout<<"ok"<<std::endl;
	}
	std::cout<<"###################"<<std::endl;
	//for(unsigned int i(0); i<yt.size(); i++){ n_box += yt[i].size();}
}

YoungTableau::~YoungTableau(){}

std::vector<YoungTableau> YoungTableau::multiply(YoungTableau const& b, unsigned int r, unsigned int c) const{
	std::vector<YoungTableau> list_tab;
	std::vector<YoungTableau> out;
	//std::cout<<"ok"<<std::endl;
	for(unsigned int k(r); k<this->row()+1; k++){
		YoungTableau tmp(*this);
		tmp.add_box(k,(1+r));
		if(tmp.valid()){ list_tab.push_back(tmp);}
	}
	//for(unsigned int k(0);k<list_tab.size();k++){
		//list_tab[k].print();
	//}
	//std::cout<<"***************"<<std::endl;
	//for(unsigned int i(0);i<list_tab.size();i++){list_tab[i].print();}
	c++;
	if(c == b.col(r)) {
		c=0;
		r++;
	}
	if(c<b.col(r) && r<b.row()) {
		//std::cout<<"again"<<std::endl;
		for(unsigned int i(0);i<list_tab.size();i++){
			std::vector<YoungTableau> tmp(list_tab[i].multiply(b,r,c));
			for(unsigned int j(0); j<tmp.size(); j++){
				bool new_tab(true);
				for(unsigned int k(0);k<out.size();k++){
					if(tmp[j] == out[k]){
						new_tab = false;
						k=out.size();
						//std::cout<<"id"<<std::endl;
					}
				}
				if(new_tab){ out.push_back(tmp[j]);}
			}
		}
		return out;
	} else { 
		return list_tab;
	}
}

void YoungTableau::add_box(unsigned int i,unsigned int j){
	if( i<this->row() ){
		yt[i].push_back(j);
	} else {
		yt.push_back(std::vector<unsigned int> (1,j));
	}
}

bool YoungTableau::valid(){
	if(this->row()!=1){
		if(boxes_in_col[0]>N){std::cout<<"p"<<std::endl;return false;}
		for(unsigned int i(1);i<boxes_in_col[0];i++){
			if(boxes_in_row[i]>boxes_in_row[i-1]){std::cout<<"q"<<std::endl;return false;}
			if(yt[i][boxes_in_row[i]-1] != 0 && 
					yt[i][boxes_in_row[i]-1] == yt[i-1][boxes_in_row[i]-1]) {std::cout<<"r"<<std::endl;return false;}
		}
	}
	return true;
}

bool YoungTableau::operator==(YoungTableau const& b){
	if(this->size() != b.size() && this->row() != b.row()){ return false; }
	for(unsigned int i(0);i<this->row();i++){
		if(this->col(i) != b.col(i)){ return false; }
		else{
			for(unsigned int j(0);j<this->col(i);j++){
				if(this->yt[i][j] != b.yt[i][j]) {return false;}
			}
		}
	}
	return true;
}

void YoungTableau::final_check(){
	if(this->row()>=N){ 
		unsigned int col_to_del(this->col(N-1));
		std::vector<std::vector<unsigned int> > tmp;
		for(unsigned int i(0); i<N-1; i++){
			unsigned int n_col(this->col(i)-col_to_del);
			if(n_col){
				tmp.push_back(std::vector<unsigned int>(n_col,0));
				for(unsigned int j(0); j<this->row()-1;j++){
					tmp[i][j] = yt[i][j+col_to_del];
					//std::cout<<k<<std::endl;
				}
			}
		}
		yt=tmp;
		for(unsigned int i(0); i<yt.size(); i++){ n_box += yt[i].size();}
	} 
}

void YoungTableau::print(){
	for(unsigned int i(0);i<boxes_in_col[0];i++){
		for(unsigned int j(0);j<boxes_in_row[i];j++){
			if(yt[i][j]){std::cout<<yt[i][j];}
			else{std::cout<<"\u2B1C";}
		}
		std::cout<<std::endl;
	}
	std::cout<<std::endl;
}

double YoungTableau::dimension(){
	double out(1.0);
	for(unsigned int i(0);i<boxes_in_col[0];i++){
		for(unsigned int j(0);j<boxes_in_row[i];j++){
			out *= (N-i+j) / ((boxes_in_row[i]-j-1.0)+(boxes_in_col[j]-i-1.0)+1.0);
		}
	}
	return out;
}

std::vector<YoungTableau> operator*(YoungTableau const& a, YoungTableau const& b){
	std::vector<YoungTableau> out(a.multiply(b));
	//for(unsigned int i(0); i<out.size(); i++){
		//out[i].final_check();
	//}
	return out;
}

