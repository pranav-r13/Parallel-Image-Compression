import numpy as np
import matplotlib.pyplot as plt
import math

# do sizes in log scale
size_peppers = 508*321
size_cameraman = 495*500
size_cookie = 2796*1414
size_reschart = 4320*2560
sz = [size_peppers, size_cameraman, size_cookie, size_reschart]
sz_log = [math.log(i) for i in sz]

# in order
sp_seq = [1.0, 1.0, 1.0, 1.0]
sp_mpi = [1.66, 1.44, 0.41, 0.19]
sp_omp = [2.24, 2.14, 2.29, 2.45]
y = [sp_seq, sp_mpi, sp_omp]

title = "Speedup vs. Image Size (log scale)"
legend = ["Sequential", "OpenMPI", "OpenMP"]

for dt in y:
    plt.plot(sz, dt)

plt.legend(legend, loc='best')
plt.xscale('log')
plt.yscale('linear')
plt.title(title)
plt.grid(True)
plt.xlabel("Image Size (Number of Pixels, Log Scale)")
plt.ylabel("Speedup (relative to sequential)")

plt.savefig('speedup_graph.png')

#plt.show()

