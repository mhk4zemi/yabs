# YABS
Yet Another Beam Solver - This toolbox helps you getting an estimation of the beam structure natural frequencies , given various environmental boundary conditions.

![image](https://github.com/user-attachments/assets/7147b24b-71e9-4a33-a7a5-40ac57d69186)



# What can you do with this ? 

&check; Get frequencies and mode shape of the 2D euler beam 

&check; Apply soil stiffness

&check; Apply Top mass and Mass moment of Inertia 

&check; Get an illustration from the beam and its selected options.

# How to run the code 


Clone the repository with the following command, make sure you have `git`installed on your system. 

```
git clone https://github.com/mhk4zemi/yabs.git
```

To run this code you will need a `Matlab` with a version `2018b` or above. 
simply run `main.m`. 

You can control all the inputs through the `Inputs.csv` file. Feel free to modify any of the paramters set in various tables in that file. 

# Planned to be released:
- [ ] Introduce geometrical softening.
- [ ] Introduce Soil P-Y curves and activate tangent/secant methods.
- [ ] preloaded frequency/modeshape calculation.
- [ ] Introduce Hydrodynamics into the mode shape calcualtion. 
- [ ] Timoshenko beam solver.


# Request features ? 
reach out to me ( mhk4zemi@gmail.com ) If you have any feature request. 

# Inputs 
Input file is a csv file , with the star table format slightly modified. The logic is re-implemented based on the repository made public in [startable-standard](https://github.com/startable/startable-standard)
