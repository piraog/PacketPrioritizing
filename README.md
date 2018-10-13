This is a simple simulation to compare different strategies to route packets with priorities in WSNs.
The algorithms implemented are described in « Scheduling for Emergency Tasks in Industrial Wireless Sensor Networks », Changqing Xia, Xi Jin, Linghe Kong and Peng Zeng 2017.

To run this simulation you will need:
	- Python 3.X
	- Numpy for your python installation

The simulation is split in 2 files:
	- Network.py that contains the class network (for graph generation)
	- Simulation.py that contains the class Simulation and import Network (for the end simulation)

To see how the networks are generated, running Network.py will display a graph and all it's attributes
You can modify the graph's parameters by modifying the call on the Network class line 201
To understand what the parameters represents, please refer to the Documentation file and the report (section 2)

To run simulations, run the Simulation.py file
It is set to plot the comparison of SFSA and OBSSA with the parameters you can see at line 227/228
You can modify the code in lines 222/265 to plot graphs or use the code lines 215/217 to produce simple results

