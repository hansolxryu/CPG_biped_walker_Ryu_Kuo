function val = get(w,prop_name)
% GET Get asset properties from the specified object
% and return the value

% This function is designed to be used with inherited classes
% by copying it to the new class and making modifications
% in the indicated section.
switch prop_name
  case 'parms' % retrieve parms from parent.parms and this.parms
    if ~isempty(w.parent)
      parentparms = get(w.(w.parent),'parms');
      pfieldnames = fieldnames(parentparms);
      pvalues = struct2cell(parentparms);
      wparms = w.parms;
      wfieldnames = fieldnames(wparms);
      wvalues = struct2cell(wparms);
      val = cell2struct( cat(1,pvalues,wvalues), cat(1,pfieldnames,wfieldnames), 1);
    else
      val = w.parms;
    end
  case 'simparms'
    if ~isempty(w.parent)
      parentparms = get(w.(w.parent),'simparms');
      pfieldnames = fieldnames(parentparms);
      pvalues = struct2cell(parentparms);
      wsimparms = w.simparms;
      wfieldnames = fieldnames(wsimparms);
      wvalues = struct2cell(wsimparms);
      val = cell2struct( cat(1,pvalues,wvalues), cat(1,pfieldnames,wfieldnames), 1);
    else
      val = w.simparms;
    end
  case 'parent'
    val = w.parent;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code in this section must be modified for inherited classes
  case 'xstar'
    val = w.xstar;
  case 'N'
    val = w.N;
  case 'gamma'
    val = w.parms.gamma;
  case 'Kp'
    val = w.parms.Kp;
  case 'P'
    val = w.parms.P;
  case 'longperiodonly'
    val = w.simparms.longperiodonly;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  otherwise
    try
      val = get(w.(w.parent), prop_name);
    catch
      error([prop_name,' Is not a valid asset property'])
    end
end
