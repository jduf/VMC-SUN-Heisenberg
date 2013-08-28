#ifndef YOUNGTABLEAU
#define YOUNGTABLEAU

#include<vector>
#include<iostream>

class YoungTableau {
	public:
		YoungTableau(std::vector<unsigned int> a);
		~YoungTableau();

		void print();
		unsigned int size() const {return n_box;}
		unsigned int row(unsigned int j=0) const {return boxes_in_row[j];}
		unsigned int col(unsigned int i) const {return boxes_in_col[i];}

		std::vector<YoungTableau> multiply(YoungTableau const& b, unsigned int r=0, unsigned int c=0) const;
		bool operator==(YoungTableau const& b);

		void final_check();
		double dimension();

	private:
		std::vector<std::vector<unsigned int> > yt;
		std::vector<unsigned int> boxes_in_row;
		std::vector<unsigned int> boxes_in_col;
		unsigned int n_box;
		unsigned int N;

		void add_box(unsigned int i,unsigned int j);
		bool valid();
};

std::vector<YoungTableau> operator*(YoungTableau const& a, YoungTableau const& b);
#endif
