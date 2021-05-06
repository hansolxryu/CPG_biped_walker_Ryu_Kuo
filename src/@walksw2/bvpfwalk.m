function xdot = bvpfwalk(t,x,unknownparms,walk,walkparms)
% BVPFWALK automatically calls fwalk in such a way that the time is normalized
% and the parameters are passed, for use by the boundary value problem
% solver, bvpsearch.
% For walksw2, there is one unknown bvp parameter, step period T:
%   xdot = bvpfwalk(t, x, T, walkparms)
% The walkparms are the structure of parameters from get(walk, 'parms').

% Art Kuo 6/2009

% For other walk models, it's possible that there will be more than one
% unknown bvp parameter to be identified by bvp4c. Also sometimes the
% unknown parameters are augmented by additional quantities, as when
% called by findgaitb or findgaitspeedb. In those cases, unknownparms
% may also contain P or Kp other other variables. So here it is necessary
% reference the first value in unknownparms.
T = unknownparms(1);
xdot = fwalk(t*T, x, walk, walkparms);
xdot = xdot*T;

end
