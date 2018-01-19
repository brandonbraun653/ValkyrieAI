%------------------------------------------------------------------------------
%
% deriv
%
% Purpose:
%   Computes the derivative of the state vector for the normalized (GM=1)
%   Kepler's problem in three dimensions
%
%------------------------------------------------------------------------------
  function [yp] = deriv (t, y)
  
  % State vector derivative
  r = y(1:3);
  v = y(4:6);
  yp = [v;-r/((norm(r))^3)];
  
  