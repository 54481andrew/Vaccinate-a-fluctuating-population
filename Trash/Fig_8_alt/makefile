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
