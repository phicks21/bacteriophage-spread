import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

class Lattice_Class:
    def __init__(self, array, lattice_length, right_most_phage, dx, dt,  rl, cross_section, phage_constant):
        """
        all units should be in terms of micrometers and seconds.
        """
        #
        self.array = array
        #the size of the lattice - array
        self.lattice_length = lattice_length
        #position of the rightmost phage (initially = 0)
        self.right_most_phage = right_most_phage
        #lattice spacing
        self.dx = dx
        #timestep
        self.dt = dt 
        #lysis rate
        self.rl = rl
        #cross section of system
        self.cross_section = cross_section
        #phage constant
        self.phage_constant = phage_constant


    def move_cells(self):
        """
        A method to move the bacteria, phage, and infected bacteria in the system.
        Overview:loop through each row in the array and sample (poisson) how many bacteria move into/out of each lattice site.
                 method requires an 'update_array' for the new values to be moved into so that cells aren't accidentally moved twice.  
        """
        update_array = self.array
        #looping through 0th row (bacteria)
        for i in range(self.lattice_length):
            #if the site population is greater than zero
            if self.array[0][i] > 0:
                #get the site population
                site_population = self.array[0][i]
                #calculate the mean number that move from random walk probability and site population
                prob_of_move = 2*((100*self.dt)/((self.dx)**2))
                #get the number that actually move from a poisson distribution centred on the mean that move
                while True:
                    #ensure this number isn't greater than the population by icking until we get a value less than the site population
                    num_that_move   = np.random.binomial(site_population,prob_of_move)
                    if num_that_move < site_population:
                        break
                
                #left most boundary condition: all move right
                if  i == 0:
                    update_array[0][i]   -= num_that_move
                    update_array[0][i+1] += num_that_move

                #right most boundary condition: all move left
                elif i == (self.lattice_length-1):
                    update_array[0][i]   -= num_that_move
                    update_array[0][i-1] += num_that_move 

                #middle ~ half move left and half move right
                else:
                    #ensure the amount we move is not greater than the lattice site
                    while True:
                        left  = np.random.binomial(num_that_move,0.5)
                        right = num_that_move - left
                        total = left + right
                        if total < site_population:
                            break

                    update_array[0][i-1] += left
                    update_array[0][i+1] += right
                    update_array[0][i]   -= total
            else:
                pass
        #do the same as above for the next two rows...
        for i in range(self.lattice_length):
            if self.array[1][i] > 0:
                site_population = int(self.array[1][i])
                prob_of_move = 2*((100*self.dt)/((self.dx)**2))
                #get the number that actually move from a poisson distribution centred on the mean that move
                while True:
                    #ensure this number isn't greater than the population by icking until we get a value less than the site population
                    num_that_move   = np.random.binomial(site_population,prob_of_move)
                    if num_that_move < site_population:
                        break
                
                if  i == 0:
                    update_array[1][i]   -= num_that_move
                    update_array[1][i+1] += num_that_move

                elif i == (self.lattice_length-1):
                    update_array[1][i]   -= num_that_move
                    update_array[1][i-1] += num_that_move 

                else:
                    while True:
                        left  = np.random.binomial(num_that_move,0.5)
                        right = num_that_move - left
                        total = left + right
                        if total < site_population:
                            break

                    update_array[1][i-1] += left
                    update_array[1][i+1] += right
                    update_array[1][i]   -= left
                    update_array[1][i]   -= right
            else:
                pass


        for i in range(self.lattice_length):
            if self.array[2][i] > 0:
                site_population = int(self.array[2][i])
                prob_of_move = 2*((100*self.dt)/((self.dx)**2))
                #get the number that actually move from a poisson distribution centred on the mean that move
                while True:
                    #ensure this number isn't greater than the population by icking until we get a value less than the site population
                    num_that_move   = np.random.binomial(site_population,prob_of_move)
                    if num_that_move < site_population:
                        break
                
                if  i == 0:
                    update_array[2][i]   -= num_that_move
                    update_array[2][i+1] += num_that_move

                elif i == (self.lattice_length-1):
                    update_array[2][i]   -= num_that_move
                    update_array[2][i-1] += num_that_move 

                else:
                    while True:
                        left  = np.random.binomial(num_that_move,0.5)
                        right = num_that_move - left
                        total = left + right
                        if total < site_population:
                            break

                    update_array[2][i-1] += left
                    update_array[2][i+1] += right
                    update_array[2][i]   -= total
            else:
                pass

        self.array = update_array

    def infect_cells(self):
        """
        a method to cause infections.
        an equation from bacteriophage ecology is used alongside poissons statistics to calculate how many move
        """
        #loop through the array
        for i in range(self.lattice_length):
            num_bacteria = self.array[0][i]
            num_phage    = self.array[1][i]
            #at each lattice site check if the number of bacteria and the number of phage are both greater than zero
            if num_phage > 0 and num_bacteria > 0:
                N = num_bacteria
                P = num_phage
                #find the mean number of infected (see bacteriophage ecology page 389 and 395)
                infection_liklihood = N*self.phage_constant*self.dt/(self.dx*self.cross_section)
                number_of_infections = np.random.binomial(P,infection_liklihood)

                #most common condition the number of infected is less than both N and P
                if number_of_infections < N and number_of_infections < P:
                    #subtract the number that infect from the population of bacteria at that site 
                    self.array[0][i] -= number_of_infections
                    #subtract the number that infect from the population of phage at that site
                    self.array[1][i] -= number_of_infections
                    #add the number that infect to the infected population at that site
                    self.array[2][i] += number_of_infections

                #less common cases - these should occur exteremely rarely - almost never for low cross sections
                #P is limiting
                elif number_of_infections < N and number_of_infections >= P:
                    #subtract the number that infect from the population of bacteria at that site 
                    self.array[0][i] -= P
                    #subtract the number that infect from the population of phage at that site
                    self.array[1][i] -= P
                    #add the number that infect to the infected population at that site
                    self.array[2][i] += P
                #N is limiting
                elif number_of_infections >= N and number_of_infections < P:
                    #subtract the number that infect from the population of bacteria at that site 
                    self.array[0][i] -= N
                    #subtract the number that infect from the population of phage at that site
                    self.array[1][i] -= N
                    #add the number that infect to the infected population at that site
                    self.array[2][i] += N
                #P and N are limiting
                elif number_of_infections >= N and number_of_infections >= P:
                    #choose one which limits most
                    if N > P:
                        #subtract the number that infect from the population of bacteria at that site 
                        self.array[0][i] -= P
                        #subtract the number that infect from the population of phage at that site
                        self.array[1][i] -= P
                        #add the number that infect to the infected population at that site
                        self.array[2][i] += P
                    elif N < P:
                        #subtract the number that infect from the population of bacteria at that site 
                        self.array[0][i] -= N
                        #subtract the number that infect from the population of phage at that site
                        self.array[1][i] -= N
                        #add the number that infect to the infected population at that site
                        self.array[2][i] += N
                    elif N == P:
                        #subtract the number that infect from the population of bacteria at that site 
                        self.array[0][i] -= N
                        #subtract the number that infect from the population of phage at that site
                        self.array[1][i] -= N
                        #add the number that infect to the infected population at that site
                        self.array[2][i] += N

            else:
                pass
    
    def burst(self):
        """
        a method for bursting infected bacteria
        """
        #loop through the third row of the array (infected row)
        for i in range(self.lattice_length):
            #get the number of infected at that ste
            num_infected = self.array[2][i]
            #if the number of infected is greater than zero
            if num_infected > 0:
                #find the mean that burst from the lysis rate times the time step times the number of infected
                prob_of_burst = self.rl*self.dt
                while True:
                    #find the number that burst ensuring that its less than the total number of infected
                    num_that_burst = np.random.binomial(num_infected, prob_of_burst)
                    if num_that_burst <= num_infected:
                        break
                #subtract the number that burst from the infected lattice site population
                self.array[2][i] -= num_that_burst
                #add the number that burst  to the corresponding phage site popuation
                self.array[1][i] += num_that_burst*100

    def find_right_most_phage(self):
        """
        a method for tracking the position of the right most phage
        it works by checking if the sum of the elements with indexes greater than the ith element is zero.
        if the sum is zero than the ith element is the location of the rightmost phage
        if the sum is greater than zero than move onto the i+1th element
        """
        for i in range(self.lattice_length):
            right_most_position = i
            right_sum = np.sum(self.array[1][i:])
            if right_sum == 0:
                break
            else:
                pass
        self.right_most_phage = right_most_position
            
