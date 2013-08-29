#ifndef YOUNGTABLEAU
#define YOUNGTABLEAU

#include<vector>
#include<iostream>

class YoungTableau {
	public:
		YoungTableau(unsigned int N);
		YoungTableau(std::vector<unsigned int> a, unsigned int N);
		~YoungTableau();

		void test();
		void print();

		std::vector<YoungTableau> multiply(YoungTableau const& b, unsigned int r=0, unsigned int c=0) const;
		bool operator==(YoungTableau const& b);

		void final_check();
		double dimension();

	private:
		std::vector<std::vector<unsigned int> > yt;
		std::vector<unsigned int> row;
		std::vector<unsigned int> col;
		unsigned int N;

		YoungTableau add_box(unsigned int row,unsigned int b_index) const;
		bool is_tableau_valid();
};

std::vector<YoungTableau> operator*(YoungTableau const& a, YoungTableau const& b);
#endif
