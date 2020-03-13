import numpy as np
import matplotlib.pyplot as plt
import sys 
import scipy

inFile1 = sys.argv[1]
inFile2 = sys.argv[2]
inFile3 = sys.argv[3]

cross_sections = []
gradients = []
errors = []

with open(inFile1, 'r') as dataFile:
    lines = dataFile.readlines()
    for line in lines:
        tokens = line.split(',')
        cross_sections.append(float(tokens[0]))


with open(inFile2, 'r') as dataFile:
    lines = dataFile.readlines()
    for line in lines:
        tokens = line.split(',')
        gradients.append(float(tokens[0]))

with open(inFile3, 'r') as dataFile:
    lines = dataFile.readlines()
    for line in lines:
        tokens = line.split(',')
        errors.append(float(tokens[0]))

plt.rcParams['figure.figsize'] = (16, 10)
plt.rcParams['font.size'] = 18
plt.plot(cross_sections,gradients)
plt.errorbar(cross_sections,gradients,errors)
plt.scatter(cross_sections,gradients)
plt.xlabel('cross_sections micrometers squared')
plt.ylabel('speed, micrometers per second')
plt.title('A Plot of the Position of the Right Most Phage in the Simulated Lattice Against Time, cross section = 50')
plt.show()

plt.plot(gradients)
plt.show()