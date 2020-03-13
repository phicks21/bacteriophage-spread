import numpy as np
import math
import matplotlib.pyplot as plt
from Lattice_Class import Lattice_Class
from matplotlib import colors

cross_sections = np.array([5,10,15,20,25,30,40,50,75,100,130,170,200,300])
gradients = []
error_on_grad = []
plt.rcParams['figure.figsize'] = (16, 10)
plt.rcParams['font.size'] = 18


for i in cross_sections:
    gradients_sub_list = []
    for k in range(10):
        times = np.loadtxt(str(i) + 'times' + str(k) + '.dat')
        positions = np.loadtxt(str(i) + 'right_most_positions' + str(k) + '.dat')
        gradient = np.polyfit(times,positions,1)[0]
        gradients_sub_list.append(gradient)
 

    gradients.append(np.average(gradients_sub_list))
    error_on_grad.append(np.std(gradients_sub_list)/np.sqrt(10))


print(gradients)
print(error_on_grad)
plt.plot(cross_sections,gradients)
plt.xlabel('cross section - micrometers^2')
plt.ylabel('speed of propogation - micrometers.s^-1')
plt.scatter(cross_sections,gradients)
plt.errorbar(cross_sections, gradients, yerr=error_on_grad)
plt.show()

np.savetxt('cross_Sections.dat', cross_sections)
np.savetxt('gradients.dat', gradients)
np.savetxt('errors_on_gradinet.dat', error_on_grad)
