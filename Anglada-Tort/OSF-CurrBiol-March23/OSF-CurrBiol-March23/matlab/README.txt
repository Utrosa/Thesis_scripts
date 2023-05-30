cReadme File
-----------

This is the Matlab code for peak extraction and simulations.


→ compute_peaks_iterations_11.m 
Is a script intending to analyze all the data of all 11 experiments except the dense rating experiment. The data is supposed to be located in  the data/ folder. Path may need to be adjusted at the top of the script. The path needs to include the folders modules (containing analysis scripts) and support-functions containing miscellaneous functions. 

→ compute_peaks_iterations_melcon_01.m
The same peak extraction for the dense rating experiment.

→ run_simulate_2models_new_01.m
Script to run the simulations for preference and interval size model in the experiment. The path may need to be adjusted at the top of the script.

→ run_simulate_2models_new_01.m
Separated script to run the combined model. The path may need to be adjusted at the top of the script.


→ do_models.m 
Run the two models. This is the main script to simulate the models and it is called the simulate_2models_new_01 script.

Interval size model:
--------------------

→ sim_aditive_model1.m
Run the interval size model, it is called by do_models.m 

Prefrence model is computed by these lines:
-------------------------------------------
  pUtil=exp(gamma*pyxi); 
  pUtil=pUtil/sum(pUtil); % normalized
  pR=pUtil';  % prior = exp (gamma*rating)
  [pRR,~]=compute_bayes_prior_model(pR,sigma_sen,p0,xxx,NK1); % analytical
  [sim_data,HAS_sim]=sim_prob_model1(data1(:,1,:),NK1,xxx,pRR,KernelWidth); % actual 

The first two lines compute the Harrison et al. model.
 compute_bayes_prior_model compute the serial reproduction model by Langlois et al. (a variant of the Griffiths and Kalish model). 
Sim_prob_model1 - is a generic function that given a conditional distribution simulates data with this conditional distribution. Since we compute all distributions analytically, we need this function to make samples.

The combined model:
------------------------
→  sim_combined_model1.m
Combined the two models in an aditive way were pereception is the input to production.