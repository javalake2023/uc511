// cppBASMastersample.cpp

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::depends(dqrng)]]

#include <Rcpp.h>
#include <RcppThread.h>
#include <dqrng.h>

using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using dqrng::dqrunif;

#include <cmath>

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


//' @name SolveCongruence
//'
//' @title Solve system of linear congruence from HIP paper to order HIP boxes.
//'
//' @description See page 5 of Robertson et al. 2018 Halton Iterative Partitioning.
//' This is essentially an internal function and shouldn't be worried about.
//'
//' @details This function was first written by Paul van Dam-Bates for the
//' package BASMasterSample.
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
//' @details This function was first written by Paul van Dam-Bates for the
//' package BASMasterSample.
//'
//' @param lxy A Matrix of lower x y coordinates of a Halton Box inside the unit box.
//' @param base Co-prime Base but generally for BAS work it is 2 or 3.
//' @param J Integer of 2 values that represent the numbers 2^J1, 3^J2.
//'
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
//' @details This function was first written by Paul van Dam-Bates for the
//' package BASMasterSample.
//'
//' @param k An integer for the starting index k >= 0.
//' @param base Co-prime Base but generally for BAS work it is 2 or 3.
//' @param n Number of samples to draw.
//'
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
//' @details This function was written by Phil Davies.
//'
//' @param boxes bla.
//' @param boxInit bla.
//' @param intB bla.
//'
//' @return A numeric vector.
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
//' @details This function was written by Phil Davies.
//'
//' @param J A numeric vector of values with which to raise the corresponding
//' element in bases to.
//' @param bases A numeric vector containing values to raised to by the corresponding
//' powers in J.
//' @param numElements The number of elements in the numeric vector bases.
//' We might dispense with this parameter at a later stage.
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
//' @details This function was first written in R by Paul van Dam-Bates for the
//' package BASMasterSample. Subsequently it was written in C/C++ by Phil Davies.
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
//' @details This function was written by Phil Davies.
//'
//' @param a Integer to find the log to base b of.
//' @param b Base
//'
//' @return The log of a to base b.
//'
//' @examples
//' # calculate log of a to base b.
//' log_a_to_base_b(2, 4)
//'
//' @export
// [[Rcpp::export(rng = false)]]
double log_a_to_base_b(long long a, int b)
{
  return log2(a) / log2(b);
}

//
//template <typename Iter, typename T>
//inline void xiota(Iter first, Iter last, T value){
//  while(first != last){
//    *first++ = value++;
//  }
//}

//template <typename T>
//inline T pop_random(std::vector<T>& v){
//  typename std::vector<T>::size_type pos = std::rand() % v.size();
//  T res = v[pos];
//  std::swap(v[pos], v.back());
//  v.pop_back();
//  return res;
//}

//// [[Rcpp::export]]
//Rcpp::IntegerVector sample_int(int n, int min, int max){
//  Rcpp::IntegerVector res(n);
//  std::vector<int> pool(max + 1 - min);
//  xiota(pool.begin(), pool.end(), min);

//  for(R_xlen_t i = 0; i < n; i++){
//    res[i] = pop_random(pool);
//  }
//  return res;
//}

//
template<typename T>
T mod(T a, int n)
{
  return a - floor(a / n) * n;
}

//void sample_int(int first, int last, std::vector<int> *out, std::size_t n, std::mt19937*g){
//  std::ranges::sample(std::views::iota(first, last), std::back_inserter(*out), n, *g);
//}



//' @name cppRSHalton
//'
//' @title Generate numbers from a Halton Sequence with a random start
//'
//' @description For efficiency, this function can generate points along a random start
//' Halton Sequence for a predefined Halton.
//'
//' @details This function was first written in R by Paul van Dam-Bates for the
//' package BASMasterSample. Subsequently it was written in C/C++ by Phil Davies.
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
  //RcppThread::Rcout << "cppRSHalton() seeds.size() : " << seeds.size() << std::endl;
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

  int b;
  int u;

  //########### Main Loop #########################################
  for (int i = 0; i < d; i++){
    b = bases[i];
    u = seeds[i];
    k1rep = rep(boxes + u, ceil(maxrep));

    int vx2end = ceil(maxrep) - 1;
    IntegerVector vx2 = Rcpp::seq(0, vx2end) * B;
    k2rep = rep_each(vx2, boxes.length());
    k = k1rep + k2rep;
    xk = mod(k, b) / b;

    for (int j = 0; j < (ceil(log_a_to_base_b(u + n, b)) + 2); j++){
      NumericVector tmp1 = floor(k / pow(b, j+1));
      int tmp000 = pow(b, (j + 2));
      NumericVector tmp2 = mod(tmp1, b) / tmp000;
      xk = xk + tmp2;
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


//
template <typename Iter, typename T>
inline void xiota(Iter first, Iter last, T value){
  while(first != last){
    *first++ = value++;
  }
}

template <typename T>
inline T pop_random(std::vector<T>& v){
  typename std::vector<T>::size_type pos = std::rand() % v.size();
  T res = v[pos];
  std::swap(v[pos], v.back());
  v.pop_back();
  return res;
}

// [[Rcpp::export]]
Rcpp::IntegerVector sample_int(int n, int min, int max){
  Rcpp::IntegerVector res(n);
  std::vector<int> pool(max + 1 - min);
  xiota(pool.begin(), pool.end(), min);

  for(R_xlen_t i = 0; i < n; i++){
    res[i] = pop_random(pool);
  }
  return res;
}


//double runif(double min, double max) {
//  return min + static_cast<double>(rand()) / RAND_MAX * (max - min);
//}


// [[Rcpp::export]]
Rcpp::NumericVector removeDuplicates(Rcpp::NumericVector vec){
  // sort vector
  std::sort(vec.begin(), vec.end());
  // remove duplicates
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
  return vec;
}


//' @name cppBASpts
//'
//' @title Generate numbers from a Halton Sequence.
//'
//' @description For efficiency, this function can generate points along a random start
//' Halton Sequence for a predefined Halton.
//'
//' @details This function was first written in R by Blair Robertson, subsequently it was
//' re-written in C/C++ by Phil Davies.
//'
//' @param n Number of points required
//' @param seeds Random starting point in each dimension
//' @param bases Co-prime base for the Halton Sequence
//'
//' @return Matrix with the columns, order of points, x in [0,1) and y in [0,1)
//'
//' @examples
//' # First 10 points in the Halton Sequence for base 2,3
//' uc511::cppBASpts(n = 10)
//' # First 10 points in the Halton Sequence for base 2,3 with
//' # starting point at the 15th and 22nd index.
//' uc511::cppBASpts(n = 10, seeds = c(14, 21))
//'
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List cppBASpts(int n = 10,
                     IntegerVector seeds = Rcpp::IntegerVector::create(),
                     NumericVector bases = Rcpp::NumericVector::create())
{
  // set default seeds values
  if (seeds.size() == 0){
    seeds = sample_int(2, 0, 62208);
  }
  // set default bases values
  if (bases.size() == 0){
    bases = {2, 3};
  }

  //RcppThread::Rcout << "cppBASpts() n     : " << n << std::endl;
  //RcppThread::Rcout << "cppBASpts() seeds : " << seeds << std::endl;
  //RcppThread::Rcout << "cppBASpts() bases : " << bases << std::endl;

  // initialise variables
  int d = bases.length();
  int u;
  int b;

  Rcpp::NumericVector xk;
  Rcpp::List          xklist;
  Rcpp::NumericMatrix pts(n, d);

  if (seeds.length() != d){
    RcppThread::Rcout << "cppBASpts() seeds.length() != d : " << std::endl;
    seeds = rep(seeds[1], d);
  }

  //########### Main Loop #########################################
  for (int i = 0; i < d; i++){
    b = bases[i];
    u = seeds[i];
    // k <- u:(u+n-1);
    Rcpp::IntegerVector ik = Rcpp::seq(u, (u+n-1));
    Rcpp::NumericVector k = as<Rcpp::NumericVector>(ik);
    xk = mod(k, b) / b;
    xklist.push_back(removeDuplicates(clone(xk)));

    for (int j = 0; j < (ceil(log_a_to_base_b(u + n, b)) + 2); j++){
      Rcpp::NumericVector tmp1 = floor(k / std::pow(b, j + 1));
      int tmp000 = std::pow(b, (j + 2));
      NumericVector tmp2 = mod(tmp1, b) / tmp000;
      xk = xk + tmp2;
    } // end for j

    // point to column i
    NumericMatrix::Column thisCol = pts(_, i);
    // propagate changes to pts.
    thisCol = thisCol + xk;

  } // end for i

  //return pts;
  return Rcpp::List::create(_["pts"]    = pts,
                            _["xklist"] = xklist);
}


//' @name cppRSHalton_br
//'
//' @title Generate numbers from a Halton Sequence with a random start
//'
//' @description For efficiency, this function can generate points along a random start
//' Halton Sequence for a predefined Halton.
//'
//' @details This function was first written in R by Paul van Dam-Bates for the
//' package BASMasterSample. Subsequently it was written in C/C++ by Phil Davies.
//'
//' @param n Number of points required
//' @param bases Co-prime base for the Halton Sequence
//' @param seeds Random starting point in each dimension
//'
//' @return Matrix with the columns, order of point, x in [0,1) and y in [0,1)
//'
//' #examples
//' # First 10 points in the Halton Sequence for base 2,3
//' # uc511::cppRSHalton_br(n = 10)
//' # First 10 points in the Halton Sequence for base 2,3 with
//' # starting point at the 15th and 22nd index.
//' # uc511::cppRSHalton_br(n = 10, seeds = c(14, 21))
//'
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List cppRSHalton_br(int n = 10,
                          NumericVector bases = Rcpp::NumericVector::create(),
                          NumericVector seeds = Rcpp::NumericVector::create())
{
  // defaults: n = 10, bases = c(2, 3), seeds = c(0, 0)

  //if (seeds.size() == 0){
  //  seeds = {0, 0};
  //}
  if (bases.size() == 0){
    bases = {2, 3};
  }

  //RcppThread::Rcout << "cppRSHalton_br() n     : " << n << std::endl;
  //RcppThread::Rcout << "cppRSHalton_br() seeds : " << seeds << std::endl;
  //RcppThread::Rcout << "cppRSHalton_br() bases : " << bases << std::endl;

  // initialize variables.
  long long UpLim = pow(10, 15);
  long long m = 0;
  long long tmpm;
  long long u2;

  int d = bases.length();
  int b;

  double min = 0.0; // Minimum value for runif
  double max = 1.0; // Maximum value for runif

  std::vector<double> u;

  Rcpp::NumericVector xk;
  Rcpp::NumericMatrix pts(n, d);
  Rcpp::NumericVector k(n);
  Rcpp::List          xklist;

  dqrng::dqRNGkind("Xoroshiro128+");
  dqrng::dqset_seed(IntegerVector::create(42));
  NumericVector xxx = dqrunif(d, min, max);
  RcppThread::Rcout << "cppRSHalton_br() xxx : " << xxx << std::endl;

  if (seeds.size() == 0){
    // seeds = numeric(d)
    seeds[d];
    // u = runif(d)
    // Generate random values using uniform distribution
    for (int i = 0; i < d; i++) {
      // Generate random uniform values (respects the current R set.seed())
      RcppThread::Rcout << "cppRSHalton_br() u : " << R::runif(min, max) << std::endl;
      u.push_back(R::runif(min, max));
    } // end for i

    // for each element of bases
    for (int i = 0; i < d; i++) {
      b = bases[i];
      m = 0;
      int j = 1;
      // while (m + (b-1)*(b^(j-1)) <= UpLim) {
      while (m + (b-1) * std::pow(b, j-1) <= UpLim){
        // m = m + (floor(u[i]*(b^j)) %% b)*(b^(j-1))
        tmpm = std::floor(u[i] * std::pow(b, j));
        m = m + (mod(tmpm, b) * std::pow(b, j-1));
        j++;
      }
      seeds[i] = m;
    } // end for i
  } // end if seeds.size()

  //########### Main Loop #########################################
  for (int i = 0; i < d; i++){
    b = bases[i];
    u2 = seeds[i];

    // Populate the vector with values from u to (u+n-1), R: k <- u:(u+n-1);
    for (int i = 0; i < n; i++) {
      k[i] = (u2 + i);
    }

    xk = mod(k, b) / b;
    xklist.push_back(removeDuplicates(clone(xk)));

    for (int j = 1; j <= (ceil(log_a_to_base_b(u2 + n, b)) + 2); j++){
      //
      NumericVector tmp1 = Rcpp::floor(k / std::pow(b, j));
      NumericVector tmp2 = mod(tmp1, b) / std::pow(b, (j + 1));
      xk = xk + tmp2;
    } // end for j

    // point to column i
    NumericMatrix::Column thisCol = pts(_, i);
    // propagate changes to pts.
    thisCol = thisCol + xk;
  } // end for i

  //return pts;
  return Rcpp::List::create(_["pts"]    = pts,
                            _["xklist"] = xklist);
}

