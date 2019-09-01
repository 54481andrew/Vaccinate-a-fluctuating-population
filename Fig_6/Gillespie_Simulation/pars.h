// Initial values for S, Ip, and P used in Sim.cpp and Sim_Base.cpp.
// All other state variables are initialized at 0. Values
// below are output from Initial/Graph_justpath.r, giving the
// values of state variables at t = 0 on the stable limit cycle. 
int SInitGlob = 429;
int IpInitGlob = 40;
int PInitGlob = 909;

// Specify parameters to be simulated over
std::vector<double> tbVals{ 120.0 };
std::vector<double> R0pVals{ 2.0 };  
std::vector<double> gampVals{ 1/(30.0) };
std::vector<double> dVals{ 1/(1.0*365.0) };
std::vector<int> npeakVals{2696};
std::vector<double> tvVals;
std::vector<double> RhoVals{ 0.1855 , 0.37}; // vaccination fraction; only used in Sim.cpp

int tvLEN = 25; // Number of vaccination times to use



