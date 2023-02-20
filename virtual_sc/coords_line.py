
import numpy as np

#make a (nr1 x nr2) array
npoints = 31

dx = np.array([1,0,0])

x0 = np.array([-15, 0,0])


output_file = open('txt_files/coords_line.txt', 'w')

for i in range(npoints):
        x= x0 + (i*dx)
        output_file.write(str(x).replace('[','').replace(']','') + '\n')


output_file.close()


