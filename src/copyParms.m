function w2 = copyParms(w1,w2)
% Copy whatever parameters inside w1 to w2 (only if w2 happens to have 
% those parameters)
    
%wTemp = w2;
allparms1 = get(w1,'parms');
allparms2 = get(w2,'parms');
parmnames1 = fieldnames(allparms1);
parmnames2 = fieldnames(allparms2);
for iparmname1 = 1:length(parmnames1)
    for iparmname2 = 1:length(parmnames2)
        if ( length(parmnames1{iparmname1}) == length(parmnames2{iparmname2}) )
            if ( parmnames1{iparmname1} == parmnames2{iparmname2} )
%                 disp('set');
%                 disp(parmnames1{iparmname1});
%                 disp(parmnames2{iparmname2});
                w2 = set(w2,parmnames1{iparmname1},get(w1,parmnames1{iparmname1}));
            end
        end
    end
end
w2 = set(w2,'xstar',get(w1,'xstar'));
end