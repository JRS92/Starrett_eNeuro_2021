// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// Generic non-inactivating conductance
// defaults reproduce a Kd channel from here:
// http://www.jneurosci.org/content/jneuro/18/7/2309.full.pdf
#ifndef HCURRENTCUSTOM
#define HCURRENTCUSTOM
#include "conductance.hpp"

//inherit conductance class spec
class HCurrentcustom: public conductance {

public:


    double m_V_half = -95.0;
    double m_V_slope = 3;

    double m_tau_A = 100;
    double m_tau_B = 1300.0;
    double m_tau_V_half = -42.2;
    double m_tau_V_slope = -8.73;


    // specify parameters + initial conditions
    HCurrentcustom(double gbar_, double E_, double m_,double m_V_half_, double m_V_slope_, double m_tau_A_, double m_tau_B_, double m_tau_V_half_, double m_tau_V_slope_)
    {
        gbar = gbar_;
        E = E_;
        m = m_;
        

        // activation parameters
        m_V_half = m_V_half_;
        m_V_slope = m_V_slope_;


        m_tau_A = m_tau_A_;
        m_tau_B = m_tau_B_;
        m_tau_V_half = m_tau_V_half_;
        m_tau_V_slope = m_tau_V_slope_;


        // defaults 
        if (isnan(gbar)) { gbar = 0; }
        
        if (isnan (E)) { E = -19.9; }


        p = 1;

        // do not allow this channel to be approximated
        
    }

    double m_inf(double, double);
    double tau_m(double, double);
    string getClass(void);
};

string HCurrentcustom::getClass(){return "HCurrentcustom";}

double HCurrentcustom::m_inf(double V, double Ca) {return 1.0/(1.0+exp((V-m_V_half)/m_V_slope));}

double HCurrentcustom::tau_m(double V, double Ca) {return m_tau_A + m_tau_B/(1+exp((V - m_tau_V_half)/m_tau_V_slope));}



#endif
