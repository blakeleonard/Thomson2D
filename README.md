# Thomson2D

A collection of five separate methods, written in MATLAB, for solving a 2D version of the [Thomson Problem.](https://en.wikipedia.org/wiki/Thomson_problem)  In this version of the problem, the goal is to determine the ground state electrostatic potential energy configuration of N number of charged particles confined to a circle of radius one.

1.) Physical Model (Euler Method)

2.) Metropolis 

3.) Simulated Annealing

4.) Replica Exchange

5.) Interacting Replicas


The evolution of complex many-body physical systems such as this are extraordinarily similar to the behavior of heuristic algorithms operating on combinatorial optimization problems.  Combinatorial optimization problems often possess a relatively large number of locally optimal pseudo-solutions, similar to the abundance of meta-stable energy states in complex physical systems.  This can make determination of the global optimum difficult, especially for heuristic algorithms which attempt to optimize a cost function locally (i.e. by iteratively making a small change in the parameters, testing the resulting change in the cost function, and allowing the change if the cost function is decreased or some other conditions are satisfied).

The concept of a local optimizer acting on some cost function in the parameter space of a combinatorial problem is mathematically equivalent to the model of a thermal system exploring its’ energy landscape.  At a certain temperature, a thermal system can realistically exchange a certain amount of energy with the environment.  This means that while the system is generally attempting to find the state of lowest energy, it can temporarily gain energy, and in so doing escape from a local energy well.  However, if the well is deeper than the realistic allowable energy gain, then the system will remain stuck in a meta-stable, locally optimal energy state indefinitely.  Analogously, if a local well of the cost function in parameter space is deeper than the realistic allowable positive gain in the cost function, then the local optimizer will remain stuck.  Note that there is a design tradeoff here.  The greater potential gain allowed in the cost function, the easier it is for the optimizer to escape potential wells, but the longer it will take to actually find a minimum because it will have a larger search space at each step, and it will be moving “uphill” more often.  This is the reason that heuristic algorithms can generally be relied upon to produce good pseudo-solutions a few percent above the optimum value, but rarely find the actual global optimum in sufficiently complex problems.

Similar to thermal systems and heuristic optimization algorithms, modeling of the physical evolution of the Thomson system from random initial conditions demonstrates that it's also susceptible to getting stuck in meta-stable energy configurations.  In fact, the 2D Thomson problem appears to be one of the simplest physical systems to exhibit this type of meta-stable behavior.  Because of  this, it was chosen as the first case for testing a novel optimization technique based on interacting replicas against some other well-established methods.  The technique was refined, applied to other problems, including the Traveling Salesman Problem, and eventually [published in the proceedings of the SAI 2016 Computing Conference.](https://arxiv.org/abs/1406.7282)

In the case of the 2D Thomson Problem, the system will begin to exhibit geometrical meta-stable states as the number of particles is increased.  At small numbers, the particles will line up along the circular boundary.  However, upon reaching a critical number, the energy can be lowered by expelling one or more particles to the center of the circle.  Ground states consisting of a certain number of particles expelled can be separated from meta-stable states consisting of a different number by very small amounts of energy, and yet the barrier between them could be enormous, as it would be the energy required to rearrange all charges.

For instance, for N = 17 particles, the ground state energy is U_0 = 133.822 with 2 expelled, while U_2 = 134.39 with 3 expelled (in units of q^2 / 4*Pi*Eps) resulting in a difference of only 0.4%.

Created by Blake Leonard

blake@leonardlabs.net





























