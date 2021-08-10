// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// GenericMH inactivating conductance
// defaults reproduce a Sodium channel from here:
// http://www.jneurosci.org/content/jneuro/18/7/2309.full.pdf
#ifndef OPTIFASTK_GHK_V2
#define OPTIFASTK_GHK_V2
#include "conductance.hpp"

//inherit conductance class spec
class optifastK_GHK_V2: public conductance {

public:


    double m_V_half = -25.5;
    double m_V_slope = -5.29;

    double h_V_half = -48.9;
    double h_V_slope = 5.18;

    double m_tau_A = 1.32;
    double m_tau_B = -1.26;
    double m_tau_V_half = -120.0;
    double m_tau_V_slope = -25.0;

    double tau_h1 = 30;
    double tau_h2 = 500;
    
    double f_1 = 0.8;
    double f_2 = 0.2;
    
    double h1 = 0;
    double h2 = 0;
    
    // specify parameters + initial conditions
    optifastK_GHK_V2(double gbar_, double E_, double m_, double h_,double g_, double h1_,double h2_, double m_V_half_, double m_V_slope_, double h_V_half_, double h_V_slope_, double m_tau_A_, double m_tau_B_, double m_tau_V_half_, double m_tau_V_slope_, double tau_h1_, double tau_h2_, double f_1_, double f_2_)
    {
        gbar = gbar_;
        E = E_;
        m = m_;
        h = h_;
        h1 = h1_;
        h2 = h2_;

        // activation parameters
        m_V_half = m_V_half_;
        m_V_slope = m_V_slope_;

        h_V_half = h_V_half_;
        h_V_slope = h_V_slope_;

        m_tau_A = m_tau_A_;
        m_tau_B = m_tau_B_;
        m_tau_V_half = m_tau_V_half_;
        m_tau_V_slope = m_tau_V_slope_;

        tau_h1 = tau_h1_;
        tau_h2 = tau_h2_;
        
        f_1 = f_1_;
        f_2 = f_2_;




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
    void integrate(double,double);
    string getClass(void);
    double getCurrent(double);
};

string optifastK_GHK_V2::getClass(){return "optifastK_GHK_V2";}

double optifastK_GHK_V2::getCurrent(double V_){
     
     if (V_ == 0){
         V_ = 0.01;}
      return g * (V_/26)* (exp((V_-E)/26) - 1)/((exp(V_/26) - 1));
 }


double optifastK_GHK_V2::m_inf(double V, double Ca) {return 1.0/(1.0+exp((V-m_V_half)/m_V_slope));}
double optifastK_GHK_V2::h_inf(double V, double Ca) {return 1.0/(1.0+exp((V-h_V_half)/h_V_slope));}
double optifastK_GHK_V2::tau_m(double V, double Ca) {return m_tau_A + m_tau_B/(1+exp((V - m_tau_V_half)/m_tau_V_slope));}

void optifastK_GHK_V2::integrate(double V, double Ca) {
    m = m + dt * (m_inf(V,Ca) - m)/tau_m(V,Ca);
    h1 =h1 + dt * (h_inf(V,Ca) - h1)/tau_h1;
    h2 =h2 + dt * (h_inf(V,Ca)-h2)/tau_h2;
    h = f_1*h1 + f_2*h2;
    g = gbar * m * h;
    
}



#endif

