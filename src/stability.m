function [eigs,vecs] = stability(w)
%  stability(w)  Performs stability analysis on a walk object w
%
%  [eigs, vecs] = stability(w) returns the eigenvalues and eigenvectors of
%    the gradient of walk object w, about its fixed point. Eigenvalues
%    are sorted in descending order by magnitude.
%  stability(w) with no output arguments displays the eigenvectors and
%    eigenvalues.

% Modified by Art (Jul 2008) to fix double display of eigenvalues 

J = gradient(w);
[vecs,dees] = eig(J);
dees = diag(dees);
[eigs,ndx] = dsort(dees);
vecs = vecs(:, ndx);


% print a report if no outputs are assigned
if nargout == 0
  %disp(eigs')
  disp(vecs)
end
