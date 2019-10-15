# CHDS
Conflict Handling in Direct Sampling


CHDS: Conflict Handling in Direct Sampling for Stochastic Simulation of Spatial Variables

This code can be used for both unconditional and conditional by DS, Conflict Handling direct sampling (CHDS), and Fast CHDS

Authors: Hesam Soltan Mohammadi, Mohammad Javad Abdollahifard, Faramarz Doulati Ardejani - 2019


Programming Language: MATLAB.


The steps to use the code:
I- Load the Training image: ti;
    The TIs used in the paper is included in this repository
II- Set the paramters:
    "params.search_radius =" : is the search radius (maximum radius of data-events)
    "params.n =" : is the maximum number of points in the data event
    "params.alpha =" : is the value of parameter 'alpha'
    "params.beta =" is the value of parameter 'beta'
    "params.m =" this parameter devides the SG to m*m windows
    "params.disp=" 0 or 1, set this parameter to one to show the progress of the simulation
    "params.simul_type=" 1, 2 or 3, this parameter defines the desired algorithm: 1- Fast CHDS 2- CHDS 3- DS
III- Form the simulation grid: 
    To do so generate a grid of 'nans' of desired size and then assign conditioning data to desired nodes.
	example:
	   Is = nan(100,150);
	   Is(18,93)=1000;
	   Is(40,33)=1010;
	The above code initializes a grid of size 100x150 with nans and assgins conditioning values of 1000 and 1010 to two particular
	locations. 
IV- call the function "do_simulation_final":
    im = do_simulation_final(Is,ti,params);
	"im" is the simulation result,




