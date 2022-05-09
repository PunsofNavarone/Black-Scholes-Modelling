#define _USE_MATH_DEFINES
/* Define _USE_MATH_DEFINES before including math.h to expose these macro
 * definitions for common math constants.  These are placed under an #ifdef
 * since these commonly-defined names are not part of the C/C++ standards.
 */

/* Definitions of useful mathematical constants
 * M_E        - e
 * M_LOG2E    - log2(e)
 * M_LOG10E   - log10(e)
 * M_LN2      - ln(2)
 * M_LN10     - ln(10)
 * M_PI       - pi
 * M_PI_2     - pi/2
 * M_PI_4     - pi/4
 * M_1_PI     - 1/pi
 * M_2_PI     - 2/pi
 * M_2_SQRTPI - 2/sqrt(pi)
 * M_SQRT2    - sqrt(2)
 * M_SQRT1_2  - 1/sqrt(2)
 

#define M_E        2.71828182845904523536
#define M_LOG2E    1.44269504088896340736
#define M_LOG10E   0.434294481903251827651
#define M_LN2      0.693147180559945309417
#define M_LN10     2.30258509299404568402
#define M_PI       3.14159265358979323846
#define M_PI_2     1.57079632679489661923
#define M_PI_4     0.785398163397448309616
#define M_1_PI     0.318309886183790671538
#define M_2_PI     0.636619772367581343076
#define M_2_SQRTPI 1.12837916709551257390
#define M_SQRT2    1.41421356237309504880
#define M_SQRT1_2  0.707106781186547524401 */

#include<iostream>
#include<math.h> /*Difference between cmath and math.h - cmath defined variables in the global namespace
and math.h in the std namespace*/


double norm_pdf (const double& x){
    return (1.0/(pow(2*M_PI, 0.5)))*exp(-0.5*x*x);
}
// This is an approximation for cdf for X~N(0,1)
// The standard normal distribution

double norm_cdf(const double& x) {
    double k = 1/(1+0.231642*x);
    double k_sum = k*(0.319381 + k*(-0.356563 + k*(1.781477 + k*(-1.821255 + 1.3302744*k))));

    if (x>= 0.0){
        return (1 - (1/(pow(2*M_PI, 0.5)))*exp(-0.5*x*x)*k_sum);
    } else {
        return 1 - norm_cdf(-x);
    }
}

double d_j(const int& j, const double& S, const double& K, const double& r,
const double& v, const double& T) {
return (log (S/K) + (r + (pow(-1 , j-1))*0.5*v*v)*T)/(v*(pow(T,0.5)));
}

// Calculate the European vanilla call price based on the parameters defined

double call_price(const double& S, const double& K, const double& r ,
const double& v , const double& T) {
return S * norm_cdf(d_j( 1 , S , K, r , v , T))-K*exp(-r*T)*
norm_cdf(d_j(2,S,K,r,v,T));
}

double put_price(const double& S, const double& K, const double& r ,
const double& v , const double& T) {
return -S * norm_cdf(-d_j( 1 , S , K, r , v , T))+K*exp(-r*T)*
norm_cdf(-d_j(2,S,K,r,v,T));
}



int main (){

double S = 1000; // Price of option
double K = 500 ; // Strike Price
double r = 0.05 ; // Interest Rate
double v = 0.25 ; // Volatility of Asset
double T = 1.0 ; // Time until expiry


double call = call_price(S , K, r , v , T) ;
double put = put_price(S , K, r , v , T) ;


std::cout << "Underlying : " << S << std::endl;
std::cout << "Strike : " << K << std::endl;
std::cout << "Risk Free Rate : " << r << std::endl;
std::cout << "Volatility : " << v << std::endl;
std::cout << "Maturity : " << T << std::endl;
std::cout << "Call Price : " <<call << std::endl;
std::cout << "Put Price : " << put << std::endl;
return 0 ;
}