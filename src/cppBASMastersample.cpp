// cppBASMastersample.cpp

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

#include <Rcpp.h>
#include <RcppThread.h>

using namespace Rcpp;

#define RCPPTHREAD_OVERRIDE_COUT 1   // std::cout override
#define RCPPTHREAD_OVERRIDE_THREAD 1 // std::thread override

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' @name cppBASMasterSample
//'
//' @title A placeholder for now.
//'
//' @description A placeholder for now.
//'
//' @returns Nothing.

//' @export
// [[Rcpp::export]]
void cppBASMasterSample() {

  RcppThread::Rcout << "cppBASMasterSample()" << std::endl;

  return;
}


//' @name SolveCongruence
//'
//' @title Solve system of linear congruence from HIP paper to order HIP boxes.
//'
//' @description See page 5 of Robertson et al. 2018 Halton Iterative Partitioning.
//' This is essentially an internal function and shouldn't be worried about.
//'
//' @param A Matrix that is in numeric for computational reasons but is the a_i solutions for all HIP boxes.
//' @param base Co-prime Base but generally for BAS work it is 2 or 3.
//' @param J Integer of 2 values that represent the numbers 2^J1, 3^J2.
//'
//' @returns A numeric vector.
//'
//' @export
// [[Rcpp::export(rng = false)]]
NumericVector SolveCongruence(NumericMatrix& A, NumericVector& base, NumericVector J)
{
  NumericVector b = NumericVector::create( pow(base[0], J[0]), pow(base[1], J[1]) );
  double B = b[0]*b[1];

  double j, ll, jj, x;
  bool test;
  NumericVector Index(A.nrow());
  double tmp;
  for(int i = 0; i < A.nrow(); i++)
  {
    j = 0; test = false;

    ll = (A(i,0) > A(i,1)) ? 0 : 1;		// for speed use the larger value to loop through.
    jj = (ll == 1) ? 0 : 1;				// Track jj as the one to test against.
    while(!test && (j < B))
    {
      x = A(i,ll) + (b[ll])*j;		// Find multiple
      tmp = floor(x / b[jj]);	// Does that particular multiple match with the other?
      test = (x - b[jj]*tmp) == A(i,jj);
      j ++;
    }
    if(test) Index(i) = x;			// Make sure there is an actual solution
    if(!test) Index(i) = -j;		// If no solution make it obvious.
  }
  return Index;
}


//' @name GetBoxIndices
//'
//' @title Fast Implementation of finding x and y order numbers
//'
//' @description Fast Implementation of finding x and y order numbers to feed into the linear congruence equation.
//' This solves equation for a_i in the HIP paper.
//' This is essentially an internal function and shouldn't be worried about.
//'
//' @param lxy A Matrix of lower x y coordinates of a Halton Box inside the unit box.
//' @param base Co-prime Base but generally for BAS work it is 2 or 3.
//' @param J Integer of 2 values that represent the numbers 2^J1, 3^J2.
//' @return A matrix of box indices
//'
//' @export
// [[Rcpp::export(rng = false)]]
NumericMatrix GetBoxIndices(NumericMatrix& lxy, IntegerVector& base, IntegerVector J)
{
  int n = lxy.nrow();
  NumericVector ai(2);
  NumericMatrix results(n, 2);

  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < 2; j++)
    {
      ai[j] = 0;
      for(int k = 1; k <= J[j]; k++)
      {
        ai[j] += (int(floor(lxy(i,j) * pow( base[j], k))) % base[j]) * pow( base[j], k - 1);
      }
      results(i,j) = ai[j];
    }
  }
  return results;
}


//' @name cppHaltonSeq
//'
//' @title Draw Halton Sequence values for a single dimension.
//'
//' @description Note that this was borrowed from the Internet and is not my implementation.
//'
//' @param k An integer for the starting index k >= 0.
//' @param base Co-prime Base but generally for BAS work it is 2 or 3.
//' @param n Number of samples to draw.
//' @return A Halton sequence of size n.
//'
//' @examples
//' uc511::cppHaltonSeq(k = 0, base = 2, n = 10)
//'
//' @export
// [[Rcpp::export(rng = false)]]
NumericVector cppHaltonSeq(const int & k, double & base, int & n)
{
  NumericVector xk(n);
  int index = k;
  double f;
  for (int i = 0; i < n; i++) {
    index = k + i;
    f = 1;
    while(index > 0){
      f = f/base;
      xk[i] = xk[i] + f*(fmod(index,base));
      index = index/base;
    }
  }
  return xk;
}


//' @name compareBoxesBoxInit
//'
//' @title title.
//'
//' @description description.
//'
//' @param boxes bla.
//' @param boxInit bla.
//' @param intB bla.
//' @return something
//'
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector compareBoxesBoxInit(Rcpp::NumericVector boxes, Rcpp::NumericVector boxInit, int intB){
  Rcpp::NumericVector newBoxes(boxes.size(), 0);
  //
  for (int i = 0; i < boxes.size(); i++){
    if (boxes[i] < boxInit[0]){
      newBoxes[i] = intB + (boxes[i] - boxInit[0]);
    } else {
      newBoxes[i] = (boxes[i] - boxInit[0]);
    }
  }
  return newBoxes;
}


//' @name cppProductPoweredElements
//'
//' @title Raise each element in a vector by a corresponding power provided in another vector,
//'        then return the product of all the results.
//'
//' @description Raise each element in a vector by a corresponding power provided in another vector,
//'        then return the product of all the results.
//'
//' @param   J           A numeric vector of values with which to raise the corresponding
//'                      element in bases to.
//' @param   bases       A numeric vector containing values to raised to by the corresponding
//'                      powers in J
//' @param   numElements The number of elements in the numeric vector bases.
//'
//' @return The product of all the powers.
//'
//' @examples
//' # calculate the product of the powered elements
//' uc511::cppProductPoweredElements(c(1, 2, 3), c(3, 2, 1), 3)
//'
//' @export
// [[Rcpp::export(rng = false)]]
int cppProductPoweredElements(NumericVector& J, NumericVector& bases, int numElements)
{
  // calculate the sum of the powered elements
  int B = 1;
  for (int i = 0; i < numElements; i++){
    int powered_value = pow(bases[i], J[i]);
    B *= powered_value;
  }
  return B;
}


//' @name cppWhere2Start
//'
//' @title Internal function to find the ordering of the first box according to the random seed.
//'
//' @description This is a function to find which Halton Box the initial BAS point from the Master Sample falls into and
//' thus use it to order the remaining boxes based on the initial. It also helps us tracks
//' the master sample index as we skip boxes that have no resource.
//'
//' @param J Definition for the number of grid cells of Halton frame.
//' @param seeds Master Sample random seed.
//' @param bases Co-prime bases should really always be 2,3
//' @param boxes ordering of boxes that have been clipped to be reordered according to the master sample seed.
//'
//' @returns vector of reordered Halton indices.
//'
//' @export
// [[Rcpp::export(rng = false)]]
NumericVector cppWhere2Start(NumericVector& J, IntegerVector& seeds, NumericVector& bases, NumericVector& boxes){

  // calculate the sum of the powered elements
  int numElements = 2; //sizeof(bases)/sizeof(bases[0]);
  int B = cppProductPoweredElements(J, bases, numElements);

  // calculate L
  numElements = seeds.length();
  IntegerVector L(numElements);
  for (int i = 0; i < numElements; i++){
    L[i] = seeds[i] % (int)(pow(bases[i], J[i]));
  }

  //boxInit <- SolveCongruence(matrix(L, ncol = 2, nrow = 1), bases, J)
  NumericMatrix mL(1, 2, L.begin());
  NumericVector boxInit = SolveCongruence(mL, bases, J);

  if(boxes.isNULL())
    return NA_REAL;

  // boxes <- ifelse(boxes < boxInit, B + (boxes - boxInit), boxes - boxInit)
  boxes = compareBoxesBoxInit(boxes, boxInit, B);
  return boxes.sort(false);
}


//' @name log_a_to_base_b
//'
//' @title Compute the log of a to base b.
//'
//' @description Compute the log of a to base b.
//'
//' @param   a       Integer to find the log to base b of.
//' @param   b       Base
//'
//' @return The log of a to base b.
//'
//' @examples
//' # calculate log of a to base b.
//' log_a_to_base_b(2, 4)
//'
//' @export
// [[Rcpp::export(rng = false)]]
double log_a_to_base_b(int a, int b)
{
  return log2(a) / log2(b);
}


//
template<typename T>
T mod(T a, int n)
{
  return a - floor(a / n) * n;
}


//' @name cppRSHalton
//'
//' @title Generate numbers from a Halton Sequence with a random start
//'
//' @description For efficiency, this function can generate points along a random start Halton Sequence for
//' predefined Halton.
//'
//' @param n Number of points required
//' @param seeds Random starting point in each dimension
//' @param bases Co-prime base for the Halton Sequence
//' @param boxes Halton boxes that points are required to be generated in
//' @param J Defines the Halton frame, and relates to the number of boxes used.
//'
//' @return Matrix with the columns, order of point, x in [0,1) and y in [0,1)
//'
//' @examples
//' # First 10 points in the Halton Sequence for base 2,3
//' uc511::cppRSHalton(n = 10)
//' # First 10 points in the Halton Sequence for base 2,3 with
//' # starting point at the 15th and 22nd index.
//' uc511::cppRSHalton(n = 10, seeds = c(14, 21))
//'
//' @export
// [[Rcpp::export(rng = false)]]
NumericVector cppRSHalton(int n = 10,
                          IntegerVector seeds = Rcpp::IntegerVector::create(),
                          NumericVector bases = Rcpp::NumericVector::create(),
                          NumericVector boxes = Rcpp::NumericVector::create(),
                          NumericVector J     = Rcpp::NumericVector::create())
{
  // defaults: n = 10, seeds = c(0,0), bases = c(2,3), boxes = 0, J = c(0,0)
  RcppThread::Rcout << "cppRSHalton() seeds.size() : " << seeds.size() << std::endl;
  if (seeds.size() == 0){
    seeds = {0, 0};
  }
  if (bases.size() == 0){
    bases = {2, 3};
  }
  if (boxes.size() == 0){
    boxes = {0};
  }
  if (J.size() == 0){
    J = {0, 0};
  }

  NumericVector xk;

  //RcppThread::Rcout << "cppRSHalton() n     : " << n << std::endl;
  //RcppThread::Rcout << "cppRSHalton() seeds : " << seeds << std::endl;
  //RcppThread::Rcout << "cppRSHalton() bases : " << bases << std::endl;
  //RcppThread::Rcout << "cppRSHalton() boxes : " << boxes << std::endl;
  //RcppThread::Rcout << "cppRSHalton() J     : " << J << std::endl;

  int d = bases.length();

  if (seeds.length() != d){
    RcppThread::Rcout << "cppRSHalton() seeds.length() != d : " << std::endl;
    seeds = rep(seeds[1], d);
  }

  NumericVector subsetJ = J[Rcpp::Range(0, 1)];
  IntegerVector subsetSeeds = seeds[Rcpp::Range(0, 1)];
  NumericVector subsetBases = bases[Rcpp::Range(0, 1)];
  boxes = cppWhere2Start(subsetJ, subsetSeeds, subsetBases, boxes);

  int B = cppProductPoweredElements(J, bases, 2);
  double maxrep = n / (double)boxes.length();
  int ceiling = ceil(maxrep);
  int repeat_count = ceiling - 1;
  int each = boxes.length();

  Rcpp::NumericVector k1rep;
  Rcpp::NumericVector k2rep;
  Rcpp::NumericVector k;

  Rcpp::NumericMatrix pts(n + ((ceiling * each) - n), 1);

  int u;

  //########### Main Loop #########################################
  for (int i = 0; i < d; i++){
    int b = bases[i];
    u = seeds[i];
    k1rep = rep(boxes + u, ceil(maxrep));

    int vx2end = ceil(maxrep) - 1;
    IntegerVector vx2 = Rcpp::seq(0, vx2end) * B;
    k2rep = rep_each(vx2, boxes.length());
    k = k1rep + k2rep;
    xk = mod(k, b) / b;

    for (int j = 0; j < (ceil(log_a_to_base_b(u + n, b)) + 2); j++){
      NumericVector tmp1 = floor(k / pow(b, j+1));
      //NumericVector tmp11 = tmp1[Rcpp::Range(0, 9)];
      //NumericVector tmp0 = mod(tmp1, b);
      //NumericVector tmp00 = tmp0[Rcpp::Range(0, 9)];
      int tmp000 = pow(b, (j + 2));
      NumericVector tmp2 = mod(tmp1, b) / tmp000;
      //NumericVector tmp3 = tmp2[Rcpp::Range(0, 9)];
      xk = xk + tmp2;
      //NumericVector tmp4 = xk[Rcpp::Range(0, 9)];
    }
    pts = Rcpp::cbind(pts, xk);
  }

  NumericVector tmpk = k + 1 - u;
  // reference the first column of pts.
  NumericMatrix::Column firstCol = pts(_, 0);
  // propagate changes to pts.
  firstCol = firstCol + tmpk;
  //pts = Rcpp::cbind(tmpk, pts);
  // need to return a matrix of (n, d+1)
  return pts;
}
