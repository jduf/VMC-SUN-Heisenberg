/*{new*/
#include "Combination.hpp"

int main(){
	Vector<unsigned int> comb;
	Combination generate;
	generate.set(2,6,comb);

	unsigned int i(0);
	do{ std::cout<<comb<<std::endl; i++; }
	while(generate.next());

	std::cout<<i<<" "<<generate.number()<<std::endl;
}
/*}*/

///*{old*/
//#include <stdio.h>
//
//void printc(int comb[], int k) {
	//printf("{");
	//int i;
	//for (i = 0; i < k; ++i)
		//printf("%d, ", comb[i]);
	//printf("\b\b}\n");
//}
//
//int next_comb(int comb[], int k, int n) {
	//int i = k - 1;
	//++comb[i];
	//while ((i >= 0) && (comb[i] >= n - k + 1 + i)) {
		//--i;
		//++comb[i];
	//}
//
	//if (comb[0] > n - k)
		//return 0;
//
	//for (i = i + 1; i < k; ++i)
		//comb[i] = comb[i - 1] + 1;
//
	//return 1;
//}
//
//int main() {
	//int n = 5;
	//int k = 3;
	//int comb[16];
//
	//int i;
	//for (i = 0; i < k; ++i)
		//comb[i] = i;
//
	//printc(comb, k);
//
	//while (next_comb(comb, k, n))
		//printc(comb, k);
//
	//return 0;
//}
///*}*/

