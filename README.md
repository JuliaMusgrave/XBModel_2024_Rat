# XBModel_2024_Rat

This git repository contains MATLAB code and data files to reproduce the figures and fitting process presented in "Analysis of metabolite and strain effects on cardiac cross-bridge dynamics using model linearisation techniques" (by Musgrave et al.). There are two main script files which can be run as is (Tested on MATLAB 2021b), provided the functions and data files are provided.

## Scripts
### plot_figures
Script which plots all of the figures from the paper. Also demonstrates the use of the final model, both in ode and linearised form.

### fitting_model_permutations
Script to facilitate the fitting of up to 64 permutations of the model to the experimental data. Provided that data are in the same format (specified in the script comments and briefly below), should be appropriate to perform the same analysis with another data set.

## Data
### rat_data.mat
Contains the experimental data used for fitting. Two variables in this mat file: 'freqs', an array containing the experimental frequencies CM was collected at; and 'data', a cell array containing all the data for the different metabolite conditions measured.
Column 1 of data contains row headings which describe the data and cols 2-6 contain the data for the 5 different metabolite conditions.

### final_fit.mat
Contains the outputs from the final fit of the model. 'Table3' and 'RMSE_improve' are tables which show the data from Table 3 and Figure 5 in the paper. 'OBJs' is a 16x4 array containing the raw objective values. 'xs' is a cell array containing the optimal parameters for each permutation of the model. 'final_params' is a vector with the optimal parameters for the final model (Table 4).

## Models
The three models included are all used typically in the plot_figures script. The [2022 model](https://doi.org/10.1016/j.mbs.2022.108922) is just included to provide plots. The linear model (XBmodel_2024_linear_perms) has been set up for multiple permutations so metabolite and strain dependence must be specified as inputs. The final ODE model (XBmodel_2024_Rat) is run with Rice_style_Fredev, but can be used for a range of different simulations - including length perturbations - with typical ode solvers
