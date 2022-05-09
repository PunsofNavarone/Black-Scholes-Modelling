#include <algorithm> 
#include <math.h>
#include <iostream>

double rand1( ) {
double x = 0.0;
double y = 0.0;
double euclid_sq = 0.0;

do {
x = 2*rand()/ static_cast<double>(RAND_MAX)-1;
y = 2*rand()/ static_cast<double>(RAND_MAX)-1;
euclid_sq = x*x + y*y;
 }
while (euclid_sq>=1) ;

return x*sqrt(-2*log(euclid_sq)/euclid_sq);
}

double monte_carlo_call_price(const int& num_sims , const double& S,
const double& K, const double& r,
const double& v, const double& T) 
{
double S_adjust= S*exp(T*(r-0.5*v*v));
double S_cur = 0;
double payoff_sum =0;

for (int i =0; i<num_sims ; i++) {
double gauss1 = rand1( );
S_cur = S_adjust * exp(sqrt(v*v*T)*gauss1);
payoff_sum += std::max(S_cur-K, 0.0);
}

return (payoff_sum/ static_cast<double>(num_sims))*exp(-r*T);
}

//Computing pricing now with Monte Carlo method
double monte_carlo_put_price(const int& num_sims, const double& S,
const double& K, const double& r,
const double& v, const double& T) {
double S_adjust = S*exp(T*(r-0.5*v*v));
double S_cur=0;
double payoff_sum = 0;


for (int i =0; i<num_sims ; i++) {
double gauss1 = rand1( );
S_cur = S_adjust * exp(sqrt(v*v*T)*gauss1);
payoff_sum += std::max(K-S_cur, 0.0);
}


return (payoff_sum/static_cast<double>(num_sims))*exp(-r*T);
}

int main () {
int num_sims = 20000000; 
double K =475;
double S =92;
double r =0.15;
double v =0.34;
double T =1;

double call=monte_carlo_call_price(num_sims, S , K, r , v , T) ;
double put=monte_carlo_put_price(num_sims, S , K, r , v , T) ;

std::cout<< "Number of Paths: " <<num_sims <<std::endl;
std::cout<< "Underlying : " << S <<std::endl;
std::cout<< "Strike: " << K <<std::endl;
std::cout<< "Risk Free Rate : " << r <<std::endl;
std::cout<< "Volatility: " << v <<std::endl;
std::cout<< "Maturity : " << T <<std::endl;
std::cout<< "Call Price: " <<call<<std::endl;
std::cout<< "Put Price: " <<put<<std::endl;

return 0;

}