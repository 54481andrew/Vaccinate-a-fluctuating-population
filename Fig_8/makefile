buildrun: Gillespie_Simulation/Sim.cpp 
	echo "Running simulations"
	g++ -o Gillespie_Simulation/main Gillespie_Simulation/Sim.cpp -fopenmp
	cd Gillespie_Simulation && time ./main
	echo
	echo
	echo "Building figure"
	Rscript Graph.r
