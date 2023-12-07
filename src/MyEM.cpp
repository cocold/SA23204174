#include <stdio.h>	
#include <Rcpp.h>
using namespace Rcpp;													
// using namespace Rcpp;
//' @title EM algorithm for the Genetic problem
//' @name MyEM 
//' @description EM algorithm for the Genetic problem
//' @param nc the initial number of c type
//' @param ni the initial value of i type
//' @param nt the initial of t type 
//' @return the estimated probability
//' @examples
//' \dontrun{
//' p=MyEM(85,96,25)
//' }
//' @export
// [[Rcpp::export]]
 NumericMatrix MyEM(int nc,int ni,int nt) 
 {
   double pc0=0.2;
   double pi0=0.3;
   double pt0=0.5;
   double cc,ci,ct;
   double ii,it,tt,T;
   double pc1,pi1,pt1,e;
   T=nc+ni+nt;
   cc=nc*pc0*pc0/(pc0*pc0+2*pc0*pi0+2*pc0*pt0);
   ci=nc*2*pc0*pi0/(pc0*pc0+2*pc0*pi0+2*pc0*pt0);
   ct=nc*2*pc0*pt0/(pc0*pc0+2*pc0*pi0+2*pc0*pt0);
   ii=ni*pi0*pi0/(2*pt0*pi0+pi0*pi0);
   it=ni*2*pi0*pt0/(2*pt0*pi0+pi0*pi0);
   tt=nt;
   pc1=(2*cc+ct+ci)/(2*T);
   pi1=(2*ii+it+ci)/(2*T);
   pt1=(2*tt+ct+it)/(2*T);
   e=(pc1-pc0)*(pc1-pc0)+(pi1-pi0)*(pi1-pi0)+(pt1-pt0)*(pt1-pt0);
   while(e>0.005)
   {
     pc0=pc1;
     pi0=pi1;
     pt0=pt1;
     cc=nc*pc0*pc0/(pc0*pc0+2*pc0*pi0+2*pc0*pt0);
     ci=nc*2*pc0*pi0/(pc0*pc0+2*pc0*pi0+2*pc0*pt0);
     ct=nc*2*pc0*pt0/(pc0*pc0+2*pc0*pi0+2*pc0*pt0);
     ii=ni*pi0*pi0/(2*pt0*pi0+pi0*pi0);
     it=ni*2*pi0*pt0/(2*pt0*pi0+pi0*pi0);
     tt=nt;
     pc1=(2*cc+ct+ci)/(2*T);
     pi1=(2*ii+it+ci)/(2*T);
     pt1=(2*tt+ct+it)/(2*T);
     e=(pc1-pc0)*(pc1-pc0)+(pi1-pi0)*(pi1-pi0)+(pt1-pt0)*(pt1-pt0);
   }
   NumericMatrix p(1,3);
   p(0,0)=pc1;
   p(0,1)=pi1;
   p(0,2)=pt1;
   return(p);
}

//' @import microbenchmark
//' @importFrom Rcpp evalCpp
//' @importFrom stats rnorm rgamma
//' @useDynLib SA23204174
