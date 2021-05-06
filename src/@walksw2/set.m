function w = set(w,varargin)
% SET Set asset properties and return the updated object
property_argin = varargin;
while length(property_argin) >= 2,
  prop = property_argin{1};
  val = property_argin{2};
  property_argin = property_argin(3:end);
  
  objparmnames = fieldnames(w.parms);       % model parameters
  if isfield(w, 'simparms')   % are there simulation parameters?
    objsimparmnames = fieldnames(w.simparms); % simulation parameters
  else
    objsimparmnames = {};
  end
  if any(strcmp(prop, objparmnames))
    w.parms.(prop) = val;
  elseif strcmp(prop, 'xstar')
    w.xstar = val;
  elseif strcmp(prop, 'N')
    w.N = val;
  elseif any(strcmp(prop, objsimparmnames))
    w.simparms.(prop) = val;
  else
    try
      w.(w.parent) = set(w.(w.parent), prop, val);
    catch
      error('Invalid %s property', class(w))
    end
  end
end
