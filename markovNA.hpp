// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// 3 State Markov Model 
// defaults reproduce a Sodium channel from here:
// http://www.jneurosci.org/content/jneuro/18/7/2309.full.pdf
#ifndef MARKOVNA
#define MARKOVNA
#include "conductance.hpp"

//inherit conductance class spec
class markovNA: public conductance {

public:

    double r2 = 0.2;
    double r4 = 0.05;
    double ra = 55;
    double sa = 6.4;
    double ka = -15.9;
    double rb = 60;
    double sb = 32;
    double kb = 10;
    double rc = 30;
    double sc = 77.5;
    double kc = 12; 
    double rd = 60;
    double sd = 32;
    double kd = 10; 
    double C = 1.0;
    double O = 0.96; 
    
    
    

    // specify parameters + initial conditions
    markovNA(double gbar_,double E_, double O_, double C_, double r2_, double r4_, double ra_, double sa_, double ka_, double rb_, double sb_, double kb_, double rc_, double sc_, double kc_, double rd_, double sd_, double kd_)
    {
        gbar = gbar_;
        O = O_;
        C = C_; 
        E = E_;

        // activation parameters
        r2 = r2_;
        r4 = r4_;
        ra = ra_;
        sa = sa_;
        ka = ka_;
        rb = rb_;
        sb = sb_;
        kb = kb_;
        rc = rc_;
        sc = sc_;
        kc = kc_;
        rd = rd_;
        sd = sd_;
        kd = kd_; 
       

        // defaults 
        if (isnan(gbar)) { gbar = 0; }
        
        
        if (isnan (E)) { E = 40; }
        
        UseMInfApproximation = 0;
        UseHInfApproximation = 0; 


        
        

        // do not allow this channel to be approximated
        
        E = 40; 
    }

    double A_inf(double, double);
    double B_inf(double, double);
    double r3(double, double);
    double r1(double, double);
    void integrate(double,double);
    double getState(int);
    string getClass(void);
 
};




string markovNA::getClass(){return "markovNA";}




double markovNA::A_inf(double V, double Ca) {return ra/(1+exp((V+sa)/ka));}
double markovNA::B_inf(double V, double Ca) {return rb/(1+exp((V+sb)/kb));}
double markovNA::r3(double V, double Ca) {return rc/(1+exp((V+sc)/kc));}
double markovNA::r1(double V, double Ca) {return rd/(1+exp((V+sd)/kd));}

void markovNA::integrate(double V, double Ca) {   
   
   O = O + dt*(A_inf(V,Ca)*C + r2*(1-C-O) - (B_inf(V,Ca) + r1(V,Ca))*O);
   C = C + dt*(B_inf(V,Ca)*O + r3(V,Ca)*(1-C-O) - (A_inf(V,Ca) + r4)*C);
   g =  gbar * fast_pow(O,3); 

}

double markovNA::getState(int idx) {
    if (idx == 1) {return O;}  //ideally this should be a switch
    else if (idx == 2) {return C;}
    else {return std::numeric_limits<double>::quiet_NaN();}
}



#endif
