function displayfields(w)
% WALKSW2/DISPLAYFIELDS Command window display of fields of a walksw2
name = w.parent;
if ~isempty(name) % display the parent if there is one
  displayfields(w.(name))
end

fprintf(1,'  xstar = [ '); fprintf(1, '%.10g ', w.xstar); fprintf(1, ']\n', w.xstar);
fprintf(1,'  gamma = %g  P = %g  Kp = %g \n', ...
      w.parms.gamma, w.parms.P, w.parms.Kp);