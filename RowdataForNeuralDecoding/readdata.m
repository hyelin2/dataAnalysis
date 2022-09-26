
%%
% spike_times=data['spike_times'] #Load spike times of all neurons
% pos=data['pos'] #Load x and y positions
% pos_times=data['pos_times'][0] #Load times at which positions were recorded

% ‘xorth’: orthonormalized posterior mean of latent variable
% 
% ‘xsm’: posterior mean of latent variable before orthonormalization
% 
% ‘Vsm’: posterior covariance between latent variables
% 
% ‘VsmGP’: posterior covariance over time for each latent variable
% 
% ‘y’: neural data used to estimate the GPFA model parameters
% "spike_times" : cell of size "number of neurons" x 1. 
%  Within spike_times{i} is a vector containing all the spike times of neuron i.
%  A continuous stream of the output variables. In this example, 
% "vels" : a matrix of size "number of recorded time points" x 2 (x and y velocities were recorded) 
%  that contains the x and y velocity components at all time points. 
% "vel_times" : is a vector that states the time at all recorded time points.