This folder contains code and datasets used in computational experiments for the paper "A branch-and-price approach for last-mile deliveries with capacitated robot stations"
by Kuzbakov Yerlan, Alfandari Laurent, Diego Delle Donne.  
To compile the code one will need to install CPLEX and to download boost library (installation directories then should be written in file CMakeLists.txt).   
One must be careful with the units in the dataset instances (coordinates and shedule time) and adjust values of truck_speed and robot_speed variables in the file global_variables.h.  
Now they are set to values assuming that coordinates are in meters and the times are in hours (Dataset I). In Dataset II coordinates are in kilometers and times are in minutes.  
The structure of the instance files:  
number of facilities, number of customers  
facilities (new line for each): id, x coord, y coord, capacity  
customers (new line for each): id, x coord, y coord, demand, delivery time
