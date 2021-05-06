function w = set(w,varargin)
% SET Set asset properties and return the updated object

% This function is designed to work well with inherited classes
% with minimal or no modification.  It must be copied to the
% appropriate class directory.

property_argin = varargin;
while length(property_argin) >= 2,
  prop = property_argin{1};
  val = property_argin{2};
  property_argin = property_argin(3:end);
  
  objparmnames = fieldnames(w.parms);
  % check to see if it's in w.parms, which means it's a property
  % of this class
  if any(strcmp(prop, objparmnames))
    w.parms.(prop) = val;
  else % or check to see if it's part of the parent's parms
    try
      w.(w.parent) = set(w.(w.parent), prop, val);
    catch
      error('Invalid %s property', class(w))
    end
  end
end
