/* 
   Initial values for S, Ip, and P used in Sim.cpp and Sim_Base.cpp.
   This state variables are initialized by an ODE system in the 
   directory Initial_Condition/, which simulates the system 
   long enough to generate stable cycles. 
*/ 

int SInitGlob = 497;
int IpInitGlob = 29;
int PInitGlob = 852;

// Specify parameters to be simulated over
std::vector<double> tbVals{ 120.0 };
std::vector<double> R0pVals{ 2.0 };  
std::vector<double> gampVals{ 1/(30.0) };
std::vector<double> dVals{ 1/(1.0*365.0) };
std::vector<int> npeakVals{2696};
std::vector<double> tvVals;
std::vector<double> NvVals{ 500.0 , 1000.0 }; // Number of vaccines, used in Sim.cpp only

int tvLEN = 25; // Number of vaccination times to use



