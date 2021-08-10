// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// persistent Sodium conductance

#ifndef NAP_V2
#define NAP_V2
#include "conductance.hpp"

//inherit conductance class spec
class NaP_V2: public conductance {

public:
    
    double m_V_half = 60;
    double m_V_slope = -5.6;
    double h_V_half = 31.6;
    double h_V_slope = 3.3;
    double h_tau_A = 67.3; 
    double h_tau_B = -27.5;
    double h_tau_C = 67.3;
    double h_tau_D = 27.5;
    double h_tau_E = 574.5;
    double h_tau_F = 62.6;

    // specify parameters + initial conditions
    NaP_V2(double gbar_, double E_, double m_, double h_, double m_V_half_, double m_V_slope_, double h_V_half_, double h_V_slope_, double h_tau_A_, double h_tau_B_, double h_tau_C_, double h_tau_D_, double h_tau_E_, double h_tau_F_)
    {
        gbar = gbar_;
        E = E_;
        m = m_;
        h = h_;
        m_V_half = m_V_half_;
        m_V_slope = m_V_slope_; 
        h_V_half = h_V_half_; 
        h_V_slope = h_V_slope_; 
        h_tau_A = h_tau_A_; 
        h_tau_B = h_tau_B_; 
        h_tau_C = h_tau_C_; 
        h_tau_D = h_tau_D_; 
        h_tau_E = h_tau_E_; 
        h_tau_F = h_tau_F_; 

        p = 1;
        q = 1;


        // defaults
        if (isnan(gbar)) { gbar = 0; }        
        if (isnan (E)) { E = 50; }
    }

    double m_inf(double, double);
    double h_inf(double, double);
    double tau_m(double, double);
    double tau_h(double, double);
    string getClass(void);
};

string NaP_V2::getClass(){return "NaP_V2";}

double NaP_V2::m_inf(double V, double Ca) {return 1.0/(1.0+exp((V+m_V_half)/m_V_slope));}
double NaP_V2::h_inf(double V, double Ca) {return 1.0/(1.0+exp((V+h_V_half)/h_V_slope));}
double NaP_V2::tau_m(double V, double Ca) {return 0.4;}
double NaP_V2::tau_h(double V, double Ca)  {return (h_tau_E/(exp((h_tau_A +V)/h_tau_B)+exp((h_tau_C+V)/h_tau_D)))+h_tau_F;}

#endif
