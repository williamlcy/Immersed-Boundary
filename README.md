# **IBFoam**
# **Description**
IBFoam is a simple **fluid-structure interaction** library that combines the **immersed boundary method** with the **Piso algorithm**.
The implementation of the immersed boundary method is based on [1].    

Since the development is not yet completed, only the work related to the immersed boundary method (2 dimensions) is uploaded. The follow-up work includes:
* Immersed Boundary Method in Three Dimensions
* Flexible module based on Euler Bernoulli beam
* Coupling of olaFlow and the immersed boundary method
* Some post-processing modules

# **Download and Compilation**
```
$ git clone https://github.com/williamlcy/Immersed-Boundary.git
```
```
$ cd solver
```
```
$ ./Allwmake
```
# **Reference**
[1] Roma, A.M., Peskin, C.S., Berger, M.J., 1999. An Adaptive Version of the Immersed Boundary Method. Journal of Computational Physics 153, 509–534. https://doi.org/10.1006/jcph.1999.6293
