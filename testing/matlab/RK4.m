%--------------------------------------------------------------------------
%
% RK4 implementation
%
%--------------------------------------------------------------------------
  function [y] = RK4(func, t, y_0, h, a, b, u)
  
  k_1 = func(t    , y_0          , a, b, u);
  k_2 = func(t+h/2, y_0+(h/2)*k_1, a, b, u);
  k_3 = func(t+h/2, y_0+(h/2)*k_2, a, b, u);
  k_4 = func(t+h  , y_0+    h*k_3, a, b, u);
  
  y = y_0 + (h/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);