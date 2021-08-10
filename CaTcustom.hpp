// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// GenericMH inactivating conductance
// defaults reproduce a Sodium channel from here:
// http://www.jneurosci.org/content/jneuro/18/7/2309.full.pdf
#ifndef CATCUSTOM
#define CATCUSTOM
#include "conductance.hpp"

//inherit conductance class spec
class CaTcustom: public conductance {

public:


    double m_V_half = -52.1;
    double m_V_slope = -5.91;

    double h_V_half = -68.1;
    double h_V_slope = 9.1;

    double m_tau_A = 26.7;
    double m_tau_B = 0;
    double m_tau_V_half = -68.1;
    double m_tau_V_slope = -20.5;


    double h_tau_A = 105.0;
    double h_tau_B = 0;
    double h_tau_V_half = 55.0;
    double h_tau_V_slope = -16.9;

    

    // specify parameters + initial conditions
    CaTcustom(double gbar_, double E_, double m_, double h_, double m_V_half_, double m_V_slope_, double h_V_half_, double h_V_slope_, double m_tau_A_, double m_tau_B_, double m_tau_V_half_, double m_tau_V_slope_, double h_tau_A_,double h_tau_B_, double h_tau_V_half_, double h_tau_V_slope_)
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

        m_tau_A = m_tau_A_;
        m_tau_B = m_tau_B_;
        m_tau_V_half = m_tau_V_half_;
        m_tau_V_slope = m_tau_V_slope_;

        h_tau_A = h_tau_A_;
        h_tau_B = h_tau_B_;
        h_tau_V_half = h_tau_V_half_;
        h_tau_V_slope = h_tau_V_slope_;

      

        // defaults 
        if (isnan(gbar)) { gbar = 0; }        
        
        if (isnan (E)) { E = 40; }

        p = 2;
        q = 1;
        
        is_calcium = true;
        // do not allow this channel to be approximated
       
    }

    double m_inf(double, double);
    double h_inf(double, double);
    double tau_m(double, double);
    double tau_h(double, double);
    string getClass(void);
   
};

string CaTcustom::getClass(){return "CaTcustom";}


double CaTcustom::m_inf(double V, double Ca) {return 1.0/(1.0+exp((V-m_V_half)/m_V_slope));}
double CaTcustom::h_inf(double V, double Ca) {return 1.0/(1.0+exp((V-h_V_half)/h_V_slope));}
double CaTcustom::tau_m(double V, double Ca) {return m_tau_A + m_tau_B/(1+exp((V - m_tau_V_half)/m_tau_V_slope));}
double CaTcustom::tau_h(double V, double Ca) {return (h_tau_A + h_tau_B/(1.0+exp((V - h_tau_V_half)/h_tau_V_slope)));}


#endif