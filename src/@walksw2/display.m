function display(wobject)
% WALKSW2/DISPLAY Command window display of a walksw2

% This function is designed so that it should work with inherited objects
% simply by copying it to the new object's directory.  The function
% displayfields will need to be modified.
for i = 1:length(wobject(:))
  w = wobject(i);
  disp(' ');
  fprintf(1,[inputname(1),' = ' class(wobject) '\n'])
  displayfields(w);
  disp(' ');
end