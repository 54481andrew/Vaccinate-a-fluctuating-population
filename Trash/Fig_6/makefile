buildrun: Gillespie_Simulation/Sim.cpp Gillespie_Simulation/Sim_Base.cpp
	echo "Running simulations without vaccination..."
	g++ -o Gillespie_Simulation/main_base Gillespie_Simulation/Sim_Base.cpp -fopenmp
	cd Gillespie_Simulation && time ./main_base
	echo
	echo
	echo "Running simulations with vaccination..."
	g++ -o Gillespie_Simulation/main Gillespie_Simulation/Sim.cpp -fopenmp
	cd Gillespie_Simulation && time ./main
	echo
	echo
	echo "Graphing..."
	Rscript Graph.r


### Example simulation times: 
# with n=50 trials
# 8: 367.25user 1.24system 0:56.66elapsed 650%CPU
# 7: 307.53user 0.04system 0:45.71elapsed 672%CPU
# 5: 255.47user 0.06system 0:59.34elapsed 430%CPU
# 4: 
# 3: 203.47user 0.00system 1:08.91
# 2: 191.20user 0.00system 1:36.24elapsed 198%CPU
# 1: 214.04user 0.01system 3:34.07elapsed 99%CPU
# 1: 263.83user 0.65system 4:24.56elapsed 99%CPU

# with n = 200 trials
#7: 1334.52user 0.18system 3:19.19elapsed 670%CPU 
#5: 1111.97user 0.35system 4:07.02elapsed 450%CPU
#1: 818.53user 0.06system 13:38.77elapsed 99%CPU 

#with n = 200 trials, verbosewritemode = true
#1232.52user 3.66system 3:33.51elapsed 578%CPU

