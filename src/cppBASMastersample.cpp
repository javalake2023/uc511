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

//' @export
// [[Rcpp::export]]
void cppBASMastersample() {

  RcppThread::Rcout << "cppBASMastersample()" << std::endl;

  return;
}

//' Solve system of linear congruence from HIP paper to order HIP boxes.
//' See page 5 of Robertson et al. 2018 Halton Iterative Partitioning.
//' This is essentially an internal function and shouldn't be worried about.
//'
//'
//' @param A Matrix that is in numeric for computational reasons but is the a_i solutions for all HIP boxes.
//' @param base Co-prime Base but generally for BAS work it is 2 or 3.
//' @param J Integer of 2 values that represent the numbers 2^J1, 3^J2.
//'
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



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//


//' Fast Implementation of finding x and y order numbers to feed into the linear congruence equation.
//' This solves equation solves for a_i in the HIP paper.
//' This is essentially an internal function and shouldn't be worried about.
//'
//'
//' @param lxy A Matrix of lower x y coordinates of a Halton Box inside the unit box.
//' @param base Co-prime Base but generally for BAS work it is 2 or 3.
//' @param J Integer of 2 values that represent the numbers 2^J1, 3^J2.
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

//' Draw Halton Sequence values for a single dimension.
//' Note that this was borrowed from the Internet and is not my implementation.
//'
//'
//' @param x An integer for the starting index k >= 0.
//' @param base Co-prime Base but generally for BAS work it is 2 or 3.
//' @param n Number of samples to draw.
//'
//' @examples
//' HaltonSeq(k = 0, base = 2, n = 10)
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



//'
//' @export
// [[Rcpp::export(rng = false)]]
NumericVector cppWhere2Start(NumericVector& J, IntegerVector& seeds, NumericVector& bases, NumericVector& boxes){

  // calculate the sum of the powered elements
  int num_elements = 2; //sizeof(bases)/sizeof(bases[0]);
  int B = 1;
  for (int i = 0; i < num_elements; i++){
    int powered_value = pow(bases[i], J[i]);
    B *= powered_value;
  }
  //NumericVector B = cumprod(NumericVector::create( pow(bases[0], J[0]), pow(bases[1], J[1]) ));
  RcppThread::Rcout << "cppWhere2Start() B : " << B << std::endl;

  //#### calc L
  num_elements = seeds.length();
  //RcppThread::Rcout << "cppWhere2Start() num_elements : " << num_elements << std::endl;
  IntegerVector L(num_elements);
  for (int i = 0; i < num_elements; i++){
    int remainder = seeds[i] % (int)(pow(bases[i], J[i]));
    //RcppThread::Rcout << "cppWhere2Start() remainder : " << remainder << std::endl;
    //const auto [q, r] = std::div(seeds[i], pow(bases[i], J[i]));
    L[i] = remainder;
    //L.push_back(remainder);
  }
  //RcppThread::Rcout << "cppWhere2Start() L : " << L << std::endl;

  //boxInit <- SolveCongruence(matrix(L, ncol = 2, nrow = 1), bases, J)
  NumericMatrix mL(1, 2, L.begin());
  NumericVector boxInit = SolveCongruence(mL, bases, J);
  //RcppThread::Rcout << "cppWhere2Start() boxInit : " << boxInit << std::endl;

  if(boxes.isNULL())
    return NA_REAL;

  // boxes <- ifelse(boxes < boxInit, B + (boxes - boxInit), boxes - boxInit)
  boxes = compareBoxesBoxInit(boxes, boxInit, B);
  //RcppThread::Rcout << "cppWhere2Start() tst boxes : " << boxes << std::endl;

  return boxes.sort(false);
}

//'
//' @export
// [[Rcpp::export(rng = false)]]
int cppSumPoweredElements(NumericVector& J, NumericVector& bases, int numElements)
{
  // calculate the sum of the powered elements
  int B = 1;
  for (int i = 0; i < numElements; i++){
    int powered_value = pow(bases[i], J[i]);
    B *= powered_value;
  }
  return B;
}


//'
//' @export
// [[Rcpp::export(rng = false)]]
NumericVector vectorMod(NumericVector k, int b){
  NumericVector xkt;
  for (int i = 0; i < k.length(); i++){
    xkt[i] = std::remainder(k[i], b);
  }
  return xkt;
}

//'
//' @export
// [[Rcpp::export(rng = false)]]
double log_a_to_base_b(int a, int b)
{
  return log2(a) / log2(b);
}

/**
 * Extend division reminder to vectors
 *
 * @param   a       Dividend
 * @param   n       Divisor
 */
template<typename T>
T mod(T a, int n)
{
  return a - floor(a / n) * n;
}



//'
//' @export
// [[Rcpp::export(rng = false)]]
NumericVector cppRSHalton(int & n, IntegerVector& seeds, NumericVector& bases, NumericVector& boxes, NumericVector& J)
{
  NumericVector xk; //(n);
  //int index = k;
  //double f;
  //for (int i = 0; i < n; i++) {
  //  index = k + i;
  //  f = 1;
  //  while(index > 0){
  //    f = f/base;
  //    xk[i] = xk[i] + f*(fmod(index,base));
  //    index = index/base;
  //  }
  //}

  RcppThread::Rcout << "cppRSHalton() n     : " << n << std::endl;
  RcppThread::Rcout << "cppRSHalton() seeds : " << seeds << std::endl;
  RcppThread::Rcout << "cppRSHalton() bases : " << bases << std::endl;
  RcppThread::Rcout << "cppRSHalton() boxes : " << boxes << std::endl;
  RcppThread::Rcout << "cppRSHalton() J     : " << J << std::endl;

  int d = bases.length();
  NumericMatrix pts(n, d);
  if (seeds.length() != d){
    RcppThread::Rcout << "cppRSHalton() seeds.length() != d : " << std::endl;
    seeds = rep(seeds[1], d);
  }

  NumericVector subsetJ = J[Rcpp::Range(0, 1)];
  IntegerVector subsetSeeds = seeds[Rcpp::Range(0, 1)];
  NumericVector subsetBases = bases[Rcpp::Range(0, 1)];
  boxes = cppWhere2Start(subsetJ, subsetSeeds, subsetBases, boxes);
  RcppThread::Rcout << "cppRSHalton() after cppWhere2Start boxes : " << boxes << std::endl;

  int B = cppSumPoweredElements(J, bases, 2);
  RcppThread::Rcout << "cppRSHalton() B : " << B << std::endl;

  //########### Just Testing ######################################

  // Second part of k.
  // ceiling(n/length(boxes)) - 1)*B
  int ceiling = ceil(n / boxes.length());
  RcppThread::Rcout << "cppRSHalton() ceiling : " << ceiling << std::endl;
  int repeat_count = ceiling - 1;
  //int step = 36; //B;
  int each = boxes.length();

  //NumericVector result(repeat_count * each);

  std::vector<int> result;
  for (int i = 0; i < ceiling; ++i) {
    int value = i * B;
    std::vector<int> temp(each, value);
    result.insert(result.end(), temp.begin(), temp.end());
  }

  // Generate the sequence using std::iota
  //std::iota(result.begin(), result.end(), 0);

  // Multiply each element by the step value
  //std::transform(result.begin(), result.end(), result.begin(),
  //               [step](int value) { return value * step; });

  Rcpp::NumericVector v2(result.begin(), result.end());
  RcppThread::Rcout << "cppRSHalton() result : " << v2.length() << std::endl;


  //Rcpp::NumericVector k1;
  Rcpp::NumericVector k1rep;
  Rcpp::NumericVector k2rep;
  Rcpp::NumericVector k;

  //########### Main Loop #########################################
  for (int i = 0; i < d; i++){
    RcppThread::Rcout << "cppRSHalton() d[i] : " << i << std::endl;
    int b = bases[i];
    int u = seeds[i];

    //k1 = boxes + u;
    k1rep = rep(boxes + u, ceil(n / boxes.length()));
    k2rep = rep(seq(0, ceil(n / boxes.length()) - 1) * B, each = boxes.length());
    k = k1rep + k2rep;
    xk = mod(k, b) / b;
    RcppThread::Rcout << "cppRSHalton() k[0] : " << k[0] << std::endl;
    RcppThread::Rcout << "cppRSHalton() b : " << b << std::endl;
    RcppThread::Rcout << "cppRSHalton() xk[0] : " << xk[0] << std::endl;

    //RcppThread::Rcout << "cppRSHalton() k1 : " << k1 << std::endl;
    RcppThread::Rcout << "cppRSHalton() k1rep.length() : " << k1rep.length() << std::endl;
    RcppThread::Rcout << "cppRSHalton() k2rep.length() : " << k2rep.length() << std::endl;

    //for (j in 1:(ceiling(logb(u+n,b)) + 2)) {
    //  xk <- xk + (floor(k/(b^j)) %% b)/(b^(j+1));
    //}
    //RcppThread::Rcout << "cppRSHalton() u : " << u << std::endl;
    //RcppThread::Rcout << "cppRSHalton() n : " << n << std::endl;
    //int jjj = ceil(log_a_to_base_b(u + n, b)) + 2;
    //RcppThread::Rcout << "cppRSHalton() jjj : " << jjj << std::endl;

    for (int j = 0; j < (ceil(log_a_to_base_b(u + n, b)) + 2); j++){
      RcppThread::Rcout << "cppRSHalton() j : " << j << std::endl;
      NumericVector tmp1 = floor(k / (b ^ j));
      NumericVector tmp2 = mod(tmp1, b);
    }
  }

  return xk;
}


