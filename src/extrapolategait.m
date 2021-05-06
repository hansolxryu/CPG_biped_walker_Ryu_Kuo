function xnew = extrapolategait(ws, w2, varargin)
% xnew = extrapolategait(w1, w2, 'parmvary', PARM) extrapolates a guess
% for a fixed point xnew, based on an array of previously known gaits w1,
% and a target gait w2 for which the fixed point is unknown.  A fit is
% performed on previous data, where PARM is the parameter of interest
% that undergoes variation.  Optional arguments include
%   'info'    Set to 1 to print out debugging info
%   'fittype' set to 1 to force linear fit, 2 for quadratic

% revised July 2008 by Art: If parmvary is not specified, extrapolategait
% tries to figure out which parameter is being varied. Also fixed a bug
% in the use of stored arrays which caused problems with parmvary1d


fittype = length(ws) - 1; % 3 values to do quadratic, 2 to do a linear fit
info = 0;
parmvary = []; % if nothing specified, we'll figure out which is being varied

opt_argin = varargin;
while length(opt_argin) >= 2,
  opt = opt_argin{1};
  val = opt_argin{2};
  opt_argin = opt_argin(3:end);
  switch opt
    case {'delta','criterion','stepsize'}
      % options for gradsearch, do nothing
    case 'info'
      info = val;
    case 'parmvary'
      parmvary = val;
    case 'fittype'
      fittype = val; % 1 for linear, 2 for quadratic
    otherwise
      error('Extrapolategait options: parmvary, info, fittype')
  end
end

% If the user didn't specify which parameter is being varied, pick the
% parameter which undergoes the largest change
if isempty(parmvary)
  parms2 = get(ws(end), 'parms'); fnames = fieldnames(parms2);
  pvalues1 = cell2mat(struct2cell(get(ws(end),'parms'))); % a numeric array of all parameter values
  pvalues2 = cell2mat(struct2cell(get(w2,'parms')));
  pdiff = abs(pvalues2 - pvalues1);
  plarge = max(abs(pvalues2),abs(pvalues1)); plarge(plarge==0) = Inf; % ignore 0 parameters
  [p, iparm] = max(pdiff ./ plarge);
  if isempty(iparm)
    error('Findgait could not determine which parameter is being varied')
  end
  parmvary = fnames{iparm};
end

% We require a minimum of two previous gaits for a linear search
% and three previous gaits for a quadratic search

% Linear extrapolation based on all previous values

%parameters = [ws.parms]
%parameters = get(ws, 'parms')
%parray = [parameters.(parmvary)]';

parray = zeros(length(ws),1);
for i = 1:length(ws)
  parray(i) = get(ws(i), parmvary);
end


% The old way: a linear prediction based on two values
% deltaparms = parray(2) - parray(1);
% newdeltaparms = w2.parms.(parmvary) - parray(2);
% deltaxs = ws(2).xstar - ws(1).xstar;
% xnew = (deltaxs/deltaparms)*newdeltaparms + ws(2).xstar

xstar = get(ws(1),'xstar');
b = zeros(length(ws), length(xstar));
for i = 1:length(ws)
  b(i,:) = get(ws(i), 'xstar');
end

if fittype < 0
  warning('extrapolategait needs at least one previous gait');
  xnew = [];
elseif fittype == 0
  xnew = get(ws,'xstar');
elseif fittype == 1
  % Linear prediction
  A = [parray ones(size(parray))];
%  b = reshape([ws.xstar], w2.N, [])';
%  xnew = [w2.parms.(parmvary) 1]*(A\b);
  xnew = [get(w2, parmvary) 1] * (A\b);
else % fittype >= 1
  % Quadratic prediction
  A = [parray.^2 parray ones(size(parray))];
  %b = reshape([ws.xstar], w2.N, [])';
  %xnew = [w2.parms.(parmvary)^2 w2.parms.(parmvary) 1]*(A\b);
  xnew = [get(w2, parmvary)^2 get(w2, parmvary) 1]*(A\b);
end

if info
  fprintf(1,'extrapolategait called with %d previous values\n', length(ws));
  fprintf(1,'  varying parameter %s\n', parmvary);
  fprintf(1,'  attempting %d order fit\n', fittype);
  A
  b
end
