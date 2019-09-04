// Initial values for S, Ip, and P used in Sim.cpp and Sim_Base.cpp.
// All other state variables are initialized at 0. Values
// below are output from Initial/Graph_justpath.r, giving the
// values of state variables at t = 0 on the stable limit cycle. 

int SInitGlob = 458;
int IpInitGlob = 253;
int PInitGlob = 0;

// Specify parameters to be simulated over
std::vector<double> tbVals{ 60.8 };
std::vector<double> R0pVals{ 1.5 };  
std::vector<double> gampVals{ 1/( 4.0*365.0) }; 
std::vector<double> dVals{ 1/( 4.0*365.0 ) };
std::vector<int> npeakVals{1107};
std::vector<double> tvVals;
std::vector<double> NvVals{ 500.0, 1000.0};
int tvLEN = 25; // Number of vaccination times to use




