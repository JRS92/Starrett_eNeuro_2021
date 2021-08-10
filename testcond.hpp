// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// GenericMH inactivating conductance
// defaults reproduce a Sodium channel from here:
// http://www.jneurosci.org/content/jneuro/18/7/2309.full.pdf
#ifndef TESTCOND
#define TESTCOND
#include "conductance.hpp"

//inherit conductance class spec
class testcond: public conductance {

public:

    

    double m_V_half = -5.54;
    double m_V_slope = -8.85;;

    double h_V_half = -39.16;
    double h_V_slope = 11.29;

    double m_tau_a = 16.39;
    double m_tau_b =  44.02;
    double m_tau_c = -0.33;
    double m_tau_d = 0.13;
    double m_tau_e = -6.77;
    double m_tau_f = 7.76; 
    double m_tau_base = 0.93; 

    double h_tau_a = 1401.65;
    double h_tau_b = 59.77;
    double h_tau_c = -6.42;
    double h_tau_d = 1.34;
    double h_tau_e = -8.26;
    double h_tau_f = 1.82; 
    double h_tau_base = 64.50; 

    // specify parameters + initial conditions
    testcond(double gbar_,double E_, double m_, double h_, double m_V_half_, double m_V_slope_, double h_V_half_, double h_V_slope_, double m_tau_a_, double m_tau_b_, double m_tau_c_, double m_tau_d_, double m_tau_e_, double m_tau_f_, double m_tau_base_, double h_tau_a_,double h_tau_b_, double h_tau_c_, double h_tau_d_, double h_tau_e_, double h_tau_f_, double h_tau_base_)
    {
        gbar = gbar_;
       
        E = E_;
        m = m_;
        h = h_;

        // activation parameters
        m_V_half = m_V_half_;
        m_V_slope = m_V_slope_;

        h_V_half = h_V_half_;
        h_V_slope = h_V_slope_;

        m_tau_a = m_tau_a_;
        m_tau_b = m_tau_b_;
        m_tau_c = m_tau_c_;
        m_tau_d = m_tau_d_;
        m_tau_e = m_tau_e_;
        m_tau_f = m_tau_f_;
        m_tau_base = m_tau_base_; 

        h_tau_a = h_tau_a_;
        h_tau_b = h_tau_b_;
        h_tau_c = h_tau_c_;
        h_tau_d = h_tau_d_;
        h_tau_e = h_tau_e_;
        h_tau_f = h_tau_f_;
        h_tau_base = h_tau_base_; 



        // defaults 
        if (isnan(gbar)) { gbar = 0; }
        
        
        if (isnan (E)) { E = -92; }


        p = 1;
        q = 1;

        // do not allow this channel to be approximated
        
        E = -92; 
    }

    double m_inf(double, double);
    double h_inf(double, double);
    double tau_m(double, double);
    double tau_h(double, double);
    string getClass(void);
     double getCurrent(double);
};




string testcond::getClass(){return "testcond";}

 double testcond::getCurrent(double V_){
     
     if (V_ == 0){
         V_ = 0.01;}
      return g * (V_/26)* (exp((V_-E)/26) - 1)/((exp(V_/26) - 1));
 }



double testcond::m_inf(double V, double Ca) {return 1.0/(1.0+exp((V-m_V_half)/m_V_slope));}
double testcond::h_inf(double V, double Ca) {return 1.0/(1.0+exp((V-h_V_half)/h_V_slope));}
double testcond::tau_m(double V, double Ca) {return m_tau_base + (m_tau_a/(1.0+exp((V+m_tau_b)/m_tau_c)))*(m_tau_d+1/(1+exp((V+m_tau_e)/m_tau_f)));}
double testcond::tau_h(double V, double Ca) {return h_tau_base + (h_tau_a/(1.0+exp((V+h_tau_b)/h_tau_c)))*(h_tau_d+1.0/(1.0+exp((V+h_tau_e)/h_tau_f)));}


#endif
