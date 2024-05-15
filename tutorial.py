## This example covers the basics of using the bode-fano module.

from circuit_reader import *

## We start by making a CircuitData object. This object is capable of creating a .cir file from a template, simulating that .cir file, and reading the results of it coming from the .txt file that is exported by wrspice.

cd = CircuitData()

##  Basic parallel circuit results
## Let's start by simulating the results from a Josephson Junction in parallel with a inductor. The parallel circuit has a second inductor that is mutually coupled with the inductor, which creates a flux in the loop. The circuit diagram (with their nodes from wrspice) looks like the following:

import matplotlib.image as mpimg  # this is just to show the circuit diagram. does not affect the rest of the code
img = mpimg.imread('examples/parallel_circuit_diagram.jpg')
imgplot = plt.imshow(img)
plt.show()

## Ground is shared between the mutually coupled circuits at node 0. There is a second, optional inductor (l1) which can be used. This value is currently set to zero, and it is possible to measure the current flowing through it without affecting the circuit itself. By using the "jj_current" variation it is also possible to see the current through the Josephson Junction using a zero-inductance inductor in parallel.

template = "parallel"

## We choose the template to use, which is "parallel" for this example. There are more templates which can be found in the templates folder, and are made to make simulations simpler by just needing to change parameter values. There are also default values coded into the template. Let's use them and see what happens.

cd.simulation_cycle(template)

## A simulation_cycle does three things: create the .cir file using default/input parameters from a template, simulate the said .cir file, and read out the results. All these can be done separately if needed. Let's look at the results from the simulation.

cd.vars
cd.measurables
cd.data

## cd.vars shows what results (measurables) we have, while cd.measurables classifies each of the measured values into phases, currents, or voltages. cd.data prints the results for all measurables through a pandas table. Let's try plotting some of the results. While we can go with matplotlib, there is also a specific plotter class that I made for this purpose.

import circuit_plotter as cplt
cplt.plot_cd_measurables(cd)

## Note that this was made to quickly look through the data, so plot limits might be not as optimal as expected.

## Now we will try making the circuit more interesting. Let's change a few parameters. First, we will change the bias current in the circuit (Give it a certain nonzero value).
cd.change_param("ibias_mag", 3e-5)

## We will also add a zero-inductance inductor next to the Josephson Junction to read its current. This is done by using a variation of the template. These variations are also currently added in a separate section of the templates.

cd.simulation_cycle(template, variation="jj_current")
cplt.plot_cd_measurables(cd)

## The currents are offset in a way such that it is hard to compare. Let's remove the offsets caused by DC bias.
cplt.plot_cd_measurables(cd, input_params="currents", remove_offsets=True)

## In this case, we notice that the Josephson Junction is in its negative inductance regime, as shown in the current being 180 degrees out of phase compared to the current in the inductor in parallel. The Josephson Junction is stable because the total inductance of the whole circuit is positive, and therefore globally stable.

## Let's change the parameters one last time. This time, we will try to induce some behavior which is no longer in the globally stable regime. Instead of changing the parameters through change_param we will import a prepared parameter file.

param_file = "params/parallel_params_example.txt"
cd.simulation_cycle(template, param_file_name=param_file, variation="jj_current")
cplt.plot_cd_measurables(cd, remove_offsets=True)

## If you'd like to look at how the weird peaks in the voltage look like, for example, you can zoom in:

cplt.plot_cd_measurables(cd, input_params="voltages", xlims=(1.030e-7, 1.033e-7))
cplt.plot_cd_measurables(cd, input_params="voltages", xlims=(1.1535e-7, 1.1565e-7))

## The voltage acts as if it were exponentially increasing, but suddenly gets damped. This is the unstable behavior that occurs when the global inductance is negative, and the system does not want to stay in that state.

## And this concludes the example! As this module is a work in progress things are subject to change, but the basic objectives of the model should stay consistent with what is shown here.

