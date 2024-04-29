# Main functions calls
I. IFmod.m is the main whole-body function call. It call upon (i) organ-specific functions: brain.m, heart.m, muscle.m, gi.m, liver.m, adipose.m; 
(ii) a function to compute fluxes which make up the right-hand side of the model ODEs, callflux.m; 
(iii) call_diet.m which select the right nutrient intake for a high-carb, high-fat, or balanced meal. The caloric intake (currently set at 800 cals) can be modified, 
as well as the breakdown of macronutrients which consititute each meal type. High-carb is (currently) 90% carbs and 10% fat; high-fat is 50% carbs and 50% fat, balanced is 70% carbs and 30% fat.

II. call_IFmod.m is the call function for IFmod.m. That's where initial conditions are set, parameter values for different modules such as exercise and diet are entered.

III. The remaining functions are routines I created for repetitive computations: 
- foldchange.m computes the fold change in a value as (new-old)/old
- driver_getdata.m sorts and organizes the output of IFmod.m per organ. It makes it easier to directly extract relevant data. Each organ is assigned a number: 1) brain, 2) heart, 3) muscle, 4) gi, 5) liver, 6) adipose, and 7) other tissues.
- Examples: Y1 is a matrix of metabolite values over time for the brain. UR3 is a matrix of uptake and release rates for each metabolite for skeletal muscle, etc.

IV. example_run.m is an example of how to run IFmod.m

# Main figure files
'figures_validation.m' reproduces the figures for the calibration and validation sections of the manuscript.
'figures_wholebody.m' reproduces the figures for the whole-body section.
'figures_organs.m' reproduces the organ-specific figures.

# HELPER_FUNCTIONS folder
Folder HELPER_FUNCTIONS contains formatting function calls for making figures. No action is required.

# SENSITIVITY_ANALYSIS folder
This is a stand-alone folder to perform the local sensitivity analysis on the model. sensitivities.m computes local sensitivity indices using a forward difference. Figures from the article can be reproduced by running 'figures_sensitivities.m'
