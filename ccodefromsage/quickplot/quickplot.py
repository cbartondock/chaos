import numpy
import matplotlib.pylab as plt
import matplotlib.cm as cm


myfile = open('jin.txt') #Open text file containing data
data = myfile.readlines() #Turn file into list of lines
myfile.close() #close file
data = [map(lambda x: float(x),line.split()) for line in data]
#The previous line does this:
#for line in data:
#    line = line.split()
#    for item in line:
#        item = float(item)

data = zip(*data) #take transpose
plt.scatter(data[0],data[4]) #scatter plot on column 0 and column 4
#these two lines set plot limits
plt.xlim(min(data[0]),max(data[0]))
plt.ylim(min(data[4]),max(data[4]))
#saves plot
plt.savefig("out.ps")
#shows plot
plt.show()


