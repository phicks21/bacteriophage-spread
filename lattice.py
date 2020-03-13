import numpy as np
import matplotlib.pyplot as plt
from Lattice_Class import Lattice_Class
from matplotlib import colors
import time


number_of_phage  = 100  #the initial number of phage
right_most_phage = 0    #the initial position of the right most phage
lattice_length   = 1000 #number of lattice sites
dx = 100 #the size of a lattice site: micrometers
dt = 10 #the timestep: seconds
rl = 0.001 #the lysis rate: events/second
phage_constant  = 0.5 # phage constant micrometers.s^-1
bacteria_conc   = 0.00001 # bacteria/micrometer^3
cross_sections   = np.array([5,10,15,20,25,30,35,40,45,50,75]) #micrometers^2
#np.array([30,50,100,200,300,600])
gradient_mean   =  []
gradient_errors = []

for i in cross_sections:
    #sweeps
    gradient_array = []
    for k in range(25):
        #number of bacteria = volume*concentration
        number_of_bacteria = int(dx*lattice_length*i*bacteria_conc)
        bacteria_array = np.zeros(lattice_length)
        #randomly assign bacteria to lattice sites
        for j in range(number_of_bacteria):
            index = np.random.randint(0,(lattice_length))
            bacteria_array[index] += 1
        #set up the matrix
        #system is (3,LatticeLength)
        #rows = lattice of bacteria (0) phage (1) and infected (3)
        #columns = lattice site populations
        phage_array = np.zeros(lattice_length)
        phage_array[0] = number_of_phage
        infected_array = np.zeros(lattice_length)
        system_array = np.vstack((bacteria_array,phage_array,infected_array))
        #create the Lattice class object
        system = Lattice_Class(system_array, lattice_length, 0, dx, dt, rl, i, phage_constant)
        print('concentrations = ' + str(i) + ' Sweep = ' + str(k))
        print(system.right_most_phage)
        right_most_positions = []
        times = []
        gradient_array = []
        #update the system 50,000 times 
        for l in range(130000):
            system.find_right_most_phage()
            right_most_positions.append(system.right_most_phage*dx)
            times.append(l*dt)
            system.move_cells()
            system.burst()
            system.infect_cells()
            #sample the right most position every 50 iterations
            if l % 1000 == 0:
                print(system.right_most_phage)
                print('iteration = ' +str(l))

            if system.right_most_phage == 999:
                break

        print(right_most_positions[-1])
        gradient = np.polyfit(times,right_most_positions,1)[0]
        gradient_array.append(gradient)

        np.savetxt(str(i) + 'times' + str(k) + '.dat', times)
        np.savetxt(str(i) + 'right_most_positions' + str(k) + '.dat', right_most_positions)


    gradient_mean.append(np.average(gradient_array))
    gradient_errors.append(np.std(gradient_array)/np.sqrt(10))
    

np.savetxt('test_gradients.dat', gradient_mean)
np.savetxt('test_cross_sections.dat', cross_sections)
np.savetxt('test_errors.dat', gradient_errors)
