
/*
This c++ script simulates a stochastic SIR epidemiological
model. This code implements the Gillespie algorithm to simulate
the use of vaccination to prevent the invasion
of a zoonotic pathogen. 
Note that in this simulation, hosts do not recover from 
pathogen infection. We use the parameter gamp to describe 
death from infection. 
*/

#include <iostream> // input/output std::cout, cin
#include <math.h>   // Math functions
#include <stdlib.h>
#include <time.h>  
#include <fstream> // Used to open files to which data are saved
#include <vector> 
#include <string.h>
#include <sys/stat.h> // Function: mkdir
#include <cmath> // Function: fmod
#include <random> // Mersenne twister RNG
#include <functional> //Used with RNG (bind)
#include "omp.h" // Allows parallel simulation
#include "pars.h" // Initial condition and parameters

char SimName[50] = "Badger";
int NTrials = 50; // Number of simulations for each parameter set

//**********************
//STORAGE ARRAYS
//**********************
//Define 2 arrays: ParMat, TExtMat. Each row of ParMat
//stores a unique set of parameters. Element i,j of
//TExtMat stores the time that the pathogen persists in
//the jth trial of parameter set i. 
int NParSets = tvLEN*tbVals.size()*R0pVals.size()*NvVals.size()*gampVals.size()*dVals.size()*npeakVals.size();
std::vector<std::vector<double>> ParMat(NParSets, std::vector<double>(11)); // NParSets x 11 matrix
std::vector<std::vector<double>> TExtMat(NParSets, std::vector<double>(NTrials)); // NParSets x NTrials matrix

//**********************
//FUNCTION DECLARATIONS
//**********************
void Initialize(); //Fills ParMat, initializes TExtMat
std::vector<double> Seq(double, double, int); //Returns sequence, used to make vector of tv values
void WriteMat(std::vector< std::vector<double> >& arr, char* filename); // Function that writes matrix to file

//Character arrays below store filenames for ParMat and TExtMat
char FileNamePar[50];
char FileNameTExt[50];
char DirName[50] = "Data/";


//**********************
//MAIN FUNCTION
//**********************
int main()
{

  //Build directory for simulation results
  strcat(DirName, SimName);  
  mkdir(DirName, ACCESSPERMS);
  strcat(DirName, "/");	  
  
  //Create filenames within directory for ParMat, TExtMat
  strcpy(FileNamePar, DirName);  
  strcat(FileNamePar, "ParMat");  
  strcpy(FileNameTExt, DirName);
  strcat(FileNameTExt, "TExtMat");
  
  Initialize(); //Fill in the parameter matrix
  WriteMat(ParMat, FileNamePar); //Write ParMat

  // pragma line splits compilation into 8 threads
#pragma omp parallel for num_threads(8)    
  for(int Par = 0; Par < NParSets; Par++) //Loop through parameter sets
    {
      
      //************************
      //USER-ASSIGNED VARIABLES
      //************************

      bool VerboseWriteFlag = true; // Write simulation data for each parameter set
      bool StopOnErad = true; // Stop simulation through trial if Ip==0
      
      int SInit = SInitGlob;
      int IpInit = IpInitGlob;
      int PInit = PInitGlob;
      
      double Nudge = 0.0000001; //Fudge parameter that moves the simulation past event conflicts (see below) 
      
      double StartTime = 0;
      double StartVaccinationTime = 2*365.0; // Vaccination occurs if t > StartVaccinationTime
      double EndTime = 8.0*365.0; //Maximum time to run each simulation trial
      double tick = 15.0; //OneSim writes data at time-intervals tick
      
      // Set up RNG via function dis(0.0,1.0): returns uniform random between 0 and 1
      std::random_device rd{};
      uint32_t seed = rd();
      int ID = omp_get_thread_num(); // Used to initialize generator
      std::mt19937 generator (ID + seed);
      std::uniform_real_distribution<double> dis(0.0, 1.0);

      // Define filenames for output from current parameter set. Initialize file with header names. 
      char FileSuffix[50], FileNameDat[50];
      strcpy(FileNameDat, DirName);
      sprintf(FileSuffix, "Par_%d",Par);
      strcat(FileNameDat, FileSuffix);
      std::ofstream out_data;
      if(VerboseWriteFlag)
	{      
	  out_data.open(FileNameDat);
	  out_data << "time S Sv Ip V P N births deaths ninfv ninfp nrecv nrecp S_death Sv_death Ip_death V_death P_death svacc npopvacc totvacc totbirthson totbirthsoff\n";
	}

      //Define simulation variables
      int S, Sv, Ip, V, P, NPop, NFails, NFailsTot, whichindex, nbirths,
	ndeaths, ninfv, ninfp, nrecv, nrecp,
	S_death, Sv_death, Ip_death, V_death, P_death, svacc, npopvacc,
	totvacc, totbirthson, totbirthsoff;
      double npeak, b0, tb, T, Nv,tv, b, gamv, R0p, Bp, gamp, d,
	Event_Rate, Event_Rate_Prod, RandDeath, dTime, t, ti,
	whichmin, unif, tmodT;
      
      //Extract parameters from ParMat
      //Parmat[0] corresponds to Par;
      b0 = ParMat[Par][1]; //b0
      d = ParMat[Par][2]; //d
      Bp = ParMat[Par][3]; //Bp
      Nv = ParMat[Par][4]; //Nv 
      tv = ParMat[Par][5]; //tv
      gamv = ParMat[Par][6]; //gamv
      gamp = ParMat[Par][7]; //gamp
      tb = ParMat[Par][8]; //tb
      T = ParMat[Par][9]; //T

      // Loop through trials for the current parameter set. 
      for(int ntrial = 0; ntrial < NTrials; ntrial++)
	{
	  //Set initial conditions
	  ti = StartTime;
	  t = StartTime;
	  S = SInit; //S
	  Sv = 0; //Sv
	  Ip = IpInit;
	  V = 0; //V
	  P = PInit; //P
	  NPop = S + Sv + Ip + V + P;
	  
	  tmodT = (double) std::fmod(t,T);
	  b = (double) b0*(tmodT < tb);
	  
	  //Track births, deaths, recoveries, etc for error checking
	  nbirths = 0; ndeaths = 0; ninfv = 0; ninfp = 0; nrecv = 0; 
	  nrecv = 0; S_death = 0; Sv_death = 0; Ip_death = 0;  
	  V_death = 0; P_death = 0; svacc = 0; npopvacc = 0; 
	  totvacc = 0; totbirthson = 0; totbirthsoff = 0;
	  
	  //Begin simulation loop. If StopOnErad==true, stop when Ip==0.
	  while(t < EndTime && (Ip > 0 || !StopOnErad))
	    {
	      {
		if( t >= ti && VerboseWriteFlag) //Write values at intervals "tick"
		  {
		    out_data << t << " " << S << " " << Sv << " " << Ip << " " << V << " " << P  << " " << NPop << " " << nbirths << " " << ndeaths <<  " " << ninfv << " " << ninfp << " " << nrecv << " " << nrecp << " " << S_death << " " << Sv_death << " " << Ip_death << " " << V_death << " " << P_death << " " << svacc << " " << npopvacc << " " << totvacc << " " << totbirthson << " " <<  totbirthsoff << "\n"; 
		    ti += tick;
		    //Reset values that are used for debugging
		    nbirths = 0; ndeaths = 0; svacc = 0; npopvacc = 0; totvacc = 0; totbirthson = 0; 
		totbirthsoff = 0; 
		  }
	  }

	      // Calculate total rate of events
	      Event_Rate = b + d*NPop + (Bp*S*Ip + Bp*Sv*Ip)/NPop + gamv*Sv + gamp*Ip;
	      Event_Rate_Prod = Event_Rate*dis(generator);
	      
	  
	      // Select random time step (doesn't allow for infite time step)
	      do
		{
		  dTime = dis(generator); //Use var dTime for storage
		  dTime = -log(dTime) / Event_Rate; //Assign next timestep
		} while (dTime == INFINITY); //Avoid infinite timesteps
	  
	  
	      // Above code gives a timestep of length dTime. However,
	      // we need to check that this does not jump over fixed-time
	      // events (start/stopping of breeding season, vaccination). 
	      // Code below calculates whichmin, the time of the next conflict
	      whichmin = T - tmodT;
	      whichindex = 0;
	      if(tmodT < tb && tb-tmodT < whichmin)
		{
		  whichmin = tb-tmodT;
		  whichindex = 1;
		}
	      if(tmodT < tv && tv-tmodT < whichmin)
		{
		  whichmin = tv-tmodT;
		  whichindex = 2;
		}

	      if(dTime < whichmin) //If no conflicts, proceed with Gillepsie event
		{
		  // 
		  if(Event_Rate_Prod <= b) //Birth
		    { S++; NPop++; nbirths++; } 
		  else if(Event_Rate_Prod <= b + d*NPop) //Death
		    {
		      RandDeath = dis(generator)*NPop;
		      if( RandDeath < S ) //S dies
			{S--; S_death++;}
		      else if(RandDeath < S + Sv) //Sv dies
			{Sv--; Sv_death++;}
		      else if(RandDeath < S + Sv + Ip) //Ip dies
			{Ip--; Ip_death++;}
		      else if(RandDeath < S + Sv + Ip + V) //V dies
			{V--; V_death++;}
		      else{P--; P_death++;}  //P dies	
		      NPop--; ndeaths++;
		    }  
		  else if(Event_Rate_Prod <= b + d*NPop + Bp*Ip*S/NPop) //Event: Pathogen infection of S
		    {S--; Ip++; ninfp++;}
		  else if(Event_Rate_Prod <= b + d*NPop + (Bp*Ip*S + Bp*Ip*Sv)/NPop ) //Event: Pathogen infection of V
		    {Sv--; Ip++; ninfp++;}
		  else if(Event_Rate_Prod <= b + d*NPop + (Bp*Ip*S + Bp*Ip*Sv)/NPop + gamv*Sv) //Event: Sv Recovery
		    {Sv--;V++;nrecv++;}
		  else if(Event_Rate_Prod <= b + d*NPop + (Bp*Ip*S + Bp*Ip*Sv)/NPop + gamv*Sv + gamp*Ip) //Event: Ip death
		    {Ip--; NPop--; } // Changed to allow death from pathogen infection. 
		}
	      else{
		//If conflict, stop at conflict and perform necessary update to state, reset
		switch(whichindex)
		  {
		  case 0 : b = b0; totbirthson++; break; //Start of birthing season
		  case 1 : b = 0.0; totbirthsoff++; break; //End of birthing season
		  case 2 : totvacc++; int temp = (int) round( (double) Nv*S/(std::max(NPop,1)));  temp  = std::min( temp,  S); temp = (t > StartVaccinationTime)*temp; S -= temp;  Sv += temp; // Annual vaccination
		  }//End else

		//Set dTime so t is updated to conflict time.
		//Nudge ensures that t moves past time of conflict
		dTime = whichmin + Nudge; 
	      } // Close if-else conflict     

	      //Update time t
	      t += dTime;
	      tmodT = std::fmod(t,T);
	    }//Close while loop 
	  
	  // Store final simulation time. If StopOnErad == true, and if the pathogen
	  // is eliminated, this is the time of elimination. 
	  TExtMat[Par][ntrial] = t;	      

	  // Write data one more time
	  if(VerboseWriteFlag){
	    out_data << t << " " << S << " " << Sv << " " << Ip << " " << V << " " << P  << " " << NPop << " " << nbirths << " " << ndeaths <<  " " << ninfv << " " << ninfp << " " << nrecv << " " << nrecp << " " << S_death << " " << Sv_death << " " << Ip_death << " " << V_death << " " << P_death << " " << svacc << " " << npopvacc << " " << totvacc << " " << totbirthson << " " << totbirthsoff << "\n"; 
	  }
	}//Close loop through NTrials
      out_data.close();
      
      //Print output
      std::cout << "*****************************" << std::endl;	  
      std::cout << "Finished Parameter Set " << Par+1 << " / " << NParSets << std::endl;
      std::cout << "*****************************" << std::endl;
    }//End Loop through NParSets (pragma omp)
  WriteMat(TExtMat, FileNameTExt); //Write TExtMat      
}//End Main







//////////////////////////////////////
///////// DEFINE FUNCTIONS ///////////
//////////////////////////////////////


//******************************************
//function to Initialize values of 2D array
//******************************************
void Initialize()
{
  double npeak, d, tb, b, gamp;
  tvVals = Seq(0.1, 364.9, tvLEN);   
  //Fill in ParMat
  int i = 0;
  for(int i1=0; i1<tvVals.size(); i1++)
    for(int i2=0; i2<R0pVals.size(); i2++) 
	  for(int i3=0; i3<NvVals.size(); i3++) 
	    for(int i4=0; i4<tbVals.size(); i4++)
	      for(int i5=0; i5<gampVals.size(); i5++)
		for(int i6=0; i6<dVals.size(); i6++)
		  for(int i7=0; i7<npeakVals.size(); i7++)
		    {
		      npeak = npeakVals[i7];
		      d = dVals[i6];
		      tb = tbVals[i4];
		      b =  d*npeak*(exp(-d*(365.0-tb))*(exp(d*365.0)-1)/(exp(d*tb)-1));
		      gamp = gampVals[i5];
		      ParMat[i][0] = i; //Par
		      ParMat[i][1] = b; //b
		      ParMat[i][2] = d; //d
		      ParMat[i][3] = R0pVals[i2]*(d + gamp); //Bp
		      ParMat[i][4] = NvVals[i3]; //Nv
		      ParMat[i][5] = tvVals[i1]; //tv
		      ParMat[i][6] = 0.07; //gamv
		      ParMat[i][7] = gampVals[i5]; //gamp
		      ParMat[i][8] = tb; //tb
		      ParMat[i][9] = 365.0; //T
		      ParMat[i][10] = npeak; //Peak population size
		      i++;
		    }
}

//************************************
//Function that returns sequence from minval to maxval, of length lengthval
//************************************
std::vector<double> Seq(double minval, double maxval, int lengthval) {
  std::vector<double> vec(lengthval);
  if(lengthval==1)
    {
      vec[0] = minval;
    }else{
  for(int seqindex = 0; seqindex < lengthval; seqindex++)
    {
      vec[seqindex] = minval + (double) seqindex/(lengthval-1)*(maxval-minval);
    }
  }
  return vec;
}


//************************************
//function to write values of 2D array
//************************************
void WriteMat(std::vector< std::vector<double> >& arr, char*filename)
{
  std::ofstream out_data;
  out_data.open(filename);
  int i,j;
  int NRow = arr.size();
  int NCol = arr[0].size();
  for(i=0; i<NRow;i++)
    {
      for(j=0; j<NCol; j++)
	{
	  	  out_data << arr[i][j] << " "; 
	}
      out_data << std::endl;
    }
  out_data.close();
}
