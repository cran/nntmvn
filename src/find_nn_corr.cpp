#include <Rcpp.h>
#include <algorithm>


using namespace Rcpp;


/*
    Find nearest neighbors based on a correlation matrix. 
    Returns an n X (m + 1) matrix. 
*/
// [[Rcpp::export]]
IntegerMatrix find_nn_corr_internal(const NumericMatrix &corrMat, int m){
	int n = corrMat.rows();
	IntegerMatrix NN(n, m + 1);
	NN.fill(NA_INTEGER);
	int *order = new int[n];
	
	for(int i = 0; i < n; i++){
		std::iota(order, order + n, 0);
		std::sort(order, order + n, [&corrMat, &i](int &j, int &k){
			return corrMat(j, i) > corrMat(k, i); });
		for(int j = 0; j < m + 1; j++)
			NN(i, j) = order[j];
		
	}
	delete[] order;
	return NN;
}
