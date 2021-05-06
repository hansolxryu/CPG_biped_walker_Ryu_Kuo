function g = gradient(w, varargin)
% J = gradient(w, varargin) Finds gradient of a walk object
% options: 'delta' (1e-4) stepsize of gradient evaluation
%          'checkxstar' (1) refine fixed point
%          'criterion' (1e-14) refined fixed point criterion

% modified by Art (7/23/2008) to revert back to original output, which was
% just the gradient instead of eigenvalues. Use stability to get that info.

delta = 1e-4;   checkxstar = 0;  criterion = 1e-13;
% If walk object has a property "ddelta" for the default delta size, e.g. because of the limitations of ADAMS, use it.
try delta = get(w,'ddelta')*10;  end;  try criterion = get(w,'dcrit');  checkxstar = 0;  end;

opt_argin = varargin;
while length(opt_argin) >= 2,
  opt = opt_argin{1};
  val = opt_argin{2};
  opt_argin = opt_argin(3:end);
  switch opt
    case 'delta'
      delta = val;
    case 'checkxstar'
      delta = val;
    case 'criterion'
      criterion = val;
    otherwise
      error('Gradient options: delta')
  end
end

if checkxstar
    w = gradsearch(w,[],'delta',delta,'criterion',criterion);
end

N = get(w,'N');

dx = delta;
dy = zeros(N,N);
ys = get(w, 'xstar'); % should just return xstar

for i=1:N,
  xn = get(w, 'xstar'); xn(i) = xn(i)+dx; dy(:,i) =  (onestep(w,xn)-ys)';
end
  
g = dy/dx;
