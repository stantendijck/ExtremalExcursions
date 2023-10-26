function storms = simulateStorms_v22(mdl, nStorms, mdlprm, data, method, simulationType, rsdModel, storm)

if nargin < 7
    rsdModel = 'none';
end
if nargin < 8 || isempty(storm)
    storm = [];
    flag = false;
else
    flag = true;
end
if length(storm) == 1
    storm = storm * ones(nStorms,1);
end
    

storms = cell(nStorms,1);
switch simulationType
    case 'peak'
        switch method
            case {'MMEM','EVAR','VAR'}
                parfor iStorm = 1:nStorms
                    if mod(iStorm,1000) == 0
                        fprintf('+');
                    elseif mod(iStorm,1000) == 0
                        fprintf('+');
                    end
                    if flag
                        storms{iStorm} = simulateStormFromPeak_v22(mdl{1}, mdl{2}, mdl{3}, mdl{4}, mdlprm, data.thr, method, rsdModel, storm(iStorm));
                    else
                        storms{iStorm} = simulateStormFromPeak_v22(mdl{1}, mdl{2}, mdl{3}, mdl{4}, mdlprm, data.thr, method, rsdModel, []);
                    end
                end
%             otherwise
%                 parfor iStorm = 1:nStorms
%                     if mod(iStorm,1000) == 0
%                         fprintf('+');
%                     elseif mod(iStorm,1000) == 0
%                         fprintf('+');
%                     end
%                     storms{iStorm} = simulateStormFromPeak_v2(mdl{1}, mdl{2}, mdl{3}, mdl{4}, mdlprm, data.thr, method, rsdModel, storm);
%                 end
        end
end
fprintf('\n');
end