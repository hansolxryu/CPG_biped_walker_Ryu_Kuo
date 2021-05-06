function val = get(w,prop_name)
% GET Get asset properties from the specified object
% and return the value

% This function is designed to be used with inherited classes
% by copying it to the new class and making minimal or no modifications
% in the indicated section.

objparmnames = fieldnames(w.parms);
if any(strcmp(prop_name, objparmnames))
  val = w.parms.(prop_name);
else
  switch prop_name
    case 'parms'
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
      case 'parent'
          val = w.parent;
      case 'testParm' %% added for test inheritance
          val = w.parms.testParm;
      case 'noiseVector'
          val = w.parms.noiseVector;
      case 'Lgain'
          val = w.parms.Lgain;
      case 'g'
          val = w.parms.g;
      case 'Mp'
          val = w.parms.Mp;
      case 'M'
          val = w.parms.M;
      case 'L'
          val = w.parms.L;
      case 'R'
          val = w.parms.R;
      case 'C'
          val = w.parms.C;
      case 'Ip'
          val = w.parms.Ip;
      case 'Il'
          val = w.parms.Il;
      case 'Tst'
          val = w.parms.Tst;
      case 'ctrlPortion'
          val = w.parms.ctrlPortion;
      case 'theta0'
          val = w.parms.theta0;
      case 'damping'
          val = w.parms.damping;
      case 'swingCtrl'
          val = w.parms.swingCtrl;
      case 'changeST'
          val = w.parms.changeST;
      case 'A0'
          val = w.parms.A0;
      case 'CutY'
          val = w.parms.CutY;
      case 'internalState' % whether GC matches
          val = w.parms.internalState;
    otherwise
      try
        val = get(w.(w.parent), prop_name);
      catch
        error([prop_name,' Is not a valid asset property'])
      end
  end
end