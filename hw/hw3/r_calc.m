function [r] = r_calc(data,v)
SS_res = sum((data - v).^2);                         % residual sum of squares
SS_tot = sum((data - mean(data)).^2 );                  % total sum of squares
r2 = 1 - (SS_res/SS_tot);                            % standard rsquared 
r = r2 * ((length(data)-1)/(length(data)-100));      % adjusted rsquared