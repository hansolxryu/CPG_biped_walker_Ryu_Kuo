function solinit = bvpsetup(walk, x0)
% BVPSETUP sets up a solution structure for boundary value solvers,
%   assuming prior existence of a limit cycle.
%
%   solinit = bvpsetup(walk, x0) takes a walk object and an initial
%     condition, integrates for one step, and puts that into the
%     solinit structure with fields x (time), y (state), and parameters
%     (a structure containing step period T).
%
%   Alternative: bvpsetup(walk) automatically uses xstar.

% Art Kuo 6/2009

if nargin == 1, % if we're not given x0, use the walk object's xstar
  x0 = get(walk, 'xstar');
end

[xc,tc,xs,ts]=onestep(walk, x0); % simulate one step and use it as a guess
solinit.x = ts' / tc; % time normalized to interval [0, 1]
solinit.y = xs';
solinit.parameters = tc; % one parameter is the guess at step period
% note that other walk classes may have more than one bvp parameter,
% and sometimes we can augment those with additional parameters such as
% when doing a findgaitspeedb search.

end