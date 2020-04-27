function [v] = v_calc(t,alpha,d,tau)
d = d * 10;
tau = tau * 10;
v = alpha * ((1/(1+exp(-t)))+exp((-t-d)/tau)-1);              
