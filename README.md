# Dissipative Flow

This program uses a finite difference scheme to find the fluid flow through a series of periodic obstacles. 

## Getting Started

### Prerequisites
* C development environment
* LAPACK; a linear algebra package for C programs (download and install here: http://www.netlib.org/lapack/)
* (Optional) MATLAB; to visualise the output flow data

To run this program, you must set up a C development environment on your machine. There are two main ways to do this:
1. Using an IDE (Integrated Development Environment), see https://fresh2refresh.com/c-programming/c-environment-setup/
2. Using a GCC compiler, see https://fresh2refresh.com/c-programming/c-environment-setup-using-gcc/

I personally prefer using the GCC method, but it is entirely up to you!

Once you have a C development environment, simply download the files from this repository.



## Running the program

You can run the program by loading the file 'flow.c' in your IDE, and then clicking 'Compile and Run'.

Alternatively you can compile it in a GCC compiler using the command

```
gcc flow.c -Wall -Werror -std=c99 -lm
```

Finally, run the MATLAB file 'plot_flow.m' to visualise the vorticity in the (*x*,*y*) plane.

![screenshot1](https://user-images.githubusercontent.com/43573338/46582470-44c03e00-ca3f-11e8-92ac-82306c6fd20e.jpg)

## How does it work?

The obstacles were assumed to be uniform in the *z* direction, so that the flow could be solved on the (*x*, *y*) plane. The obstacles also provide a frictional force that slows the fluid. Overall, the average fluid flow is constant in time.

The equations that had to be solved were the incompressible Navier-Stokes equations fo the fluid velovity **u**, which only has components in the (*x*, *y*) plane. These are solved by time-evolving the scalar fluid vorticity ![](http://latex.codecogs.com/gif.latex?%5Comega%3D%28%5Cnabla%5Ctimes%20%5Cbold%20u%29%5Ccdot%20%5Chat%7B%5Cbold%20z%7D) as a function of space and time:

![](http://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cpartial%20%5Comega%7D%7B%5Cpartial%20t%7D&plus;%5Cbold%20u%5Ccdot%20%5Cnabla%5Comega%20%3D%20%5Chat%7B%5Cbold%20z%7D%5Ccdot%20%28%5Cnabla%5Ctimes%5Cbold%20F%29&plus;%5Cnu%20%5Cnabla%5E2%5Comega)

where **F**(*x*, *y*, *t*) is the external force, and *ν* is the viscosity term.

To account for the obstacles and a constant external driving force, the external force is

![](http://latex.codecogs.com/gif.latex?%5Cbold%20F%3D%5Cbold%20G-K%28x%2Cy%29%5Cbold%20u)

where **G** is the time-varying constant, and **K**(*x*, *y*) is the frictional force due to the obstacles.

To find the velocity **u** from *ω*, a fluid potential field *φ* is used. Then, the velocity can be written as ![](http://latex.codecogs.com/gif.latex?%5Cbo%5Cbold%20u%20%3D%20%28%5Cnabla%20%5Cphi%29%5Ctimes%20%5Chat%7B%5Cbold%20z%7D&plus;%5Cbold%20u_0), and in this program, **u**<sub>0</sub> is set by the input parameters. To find *φ*, the following equation is solved:

![](http://latex.codecogs.com/gif.latex?%5Chat%7B%5Cbold%20z%7D%5Ccdot%20%28%5Cnabla%5Ctimes%5Cbold%20u%29%3D-%5Cnabla%5E2%5Cphi%3D%5Comega)

The boundary conditions in the *x* and *y* direction are periodic, with

![](http://latex.codecogs.com/gif.latex?%5Cbold%20u%28x%2Cy%29%3D%5Cbold%20u%28x&plus;IL_x%2Cy&plus;JL_y%29)

for integers *I* and *J*. Furthermore, the vorticity at *t*=0 is assumed to be zero everywhere.

So, given a vorticity field *ω*(*x*, *y*), a solution scheme requires iterating three steps:

1. Calculate the flow field *φ* given the scalar vorticity *ω*
2. Calculate the velocity **u**
3. Use this velocity to calculate the vorticity field at the next timestep.


The program was written with reference to the following specification:

**Specification:**
The input and output *x* and *y* grid have *M* and *N* grid points respectively, including the end of the domain:

*x*<sub>0</sub>=0 to *x*<sub>*M*-1</sub>=L<sub>*x*</sub>-*δx*

and

*y*=0 to *y*<sub>*N*-1</sub>=*L*<sub>*y*</sub>-*δy*.

**Input:** The input parameters are read from a file called 'input.txt', containing
1. *L<sub>x</sub>*: right *x* boundary of domain
2. *L<sub>y</sub>*: right *y* boundary of domain
3. *M*: number of *x* grid points
4. *N*: number of *y* grid points
5. *t<sub>f</sub>* : time to run simulation over
6. *t<sub>d</sub>*: simulation diagnostic timestep
7. *ν*: viscosity
8. *u*<sub>0*x*</sub>: *x* component of average velocity
9. *u*<sub>0*x*</sub>: *y* component of average velocity

The code outputs the vorticity *ω*(*x*, *y*) to a file called ‘output.txt’ as a single column, for all the grid points (*x*<sub>*i*</sub>, *y*<sub>*j*</sub>) for each time.

## Acknowledgments
* This program was for an academic assignment created by Prof. B. McMillan, Department of Physics, University of Warwick.
* The code used for manipulating 2-dimensional matrices was also acquired from Prof. B. McMillan, Department of Physics, University of Warwick.
