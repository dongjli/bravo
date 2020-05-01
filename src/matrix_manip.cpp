#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector colSumSq_dge(NumericVector x,IntegerVector dim) {
  int n = dim[0];
  int p = dim[1];
  NumericVector y(p);
  NumericVector::iterator it = x.begin();
  for(int j=0; j<p;++j)
  {
    double sumsq = 0.0;
    for(int i=0;i<n;++i)
    {
      sumsq += (*it)*(*it);
      ++it;
    }
    y[j] = sumsq;
  }
  return(y);
}

// [[Rcpp::export]]
NumericVector colSumSq_matrix(NumericMatrix x) {
  int n = x.nrow();
  int p = x.ncol();
  NumericVector y(p);
  NumericVector::iterator it = x.begin();
  for(int j=0; j<p;++j)
  {
    double sumsq = 0.0;
    for(int i=0;i<n;++i)
    {
      sumsq += (*it)*(*it);
      ++it;
    }
    y[j] = sumsq;
  }
  return(y);
}

// [[Rcpp::export]]
NumericVector colMSD_dgc(S4 mat,NumericVector m) {
  IntegerVector dims = mat.slot("Dim");
  int ncol = dims[1];
  int nrow = dims[0];

  int *p = INTEGER(mat.slot("p"));
  double *x = REAL(mat.slot("x"));

  NumericVector y(ncol);
  for(int j=0;j<ncol;++j)
  {
    double ssq = 0.0;
    double mj = m[j];
    int start = p[j];
    int end = p[j+1];
    for(int i=start;i<end;++i)
    {
      double temp = x[i] - mj;
      ssq += temp*temp;
    }
    y[j] = (ssq + mj*mj*(nrow - (end - start + 0.0)))/(nrow-1);
  }
 return(y);
}

// [[Rcpp::export]]
NumericVector colSUMIDX_dgc(S4 mat) {
  IntegerVector dims = mat.slot("Dim");
  int ncol = dims[1];

  int *p = INTEGER(mat.slot("p"));
  int *idx = INTEGER(mat.slot("i"));
  double *x = REAL(mat.slot("x"));

  NumericVector y(ncol);
  for(int j=0;j<ncol;++j)
  {
    int rowidx_sum = 0;
    int start = p[j];
    int end = p[j+1];
    for(int i=start;i<end;++i)
    {
      if(x[i] == 1) {
        rowidx_sum += idx[i];
      }	
    }
    y[j] =  rowidx_sum;
  }
 return(y);
}
