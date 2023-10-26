function storms = mysim2(mdl, nStorms, methods, simulationType, rsdModel, data, peaks)

maxk = 4; maxp = 4;
% maxk = 10; maxp = 10;
storms = struct();
for iMethod = 1:length(methods)
    method = methods{iMethod};
    switch method
        case {'MEM','MMEM'}
%             maxk = 6;
            p = 1;
%             storms.(method) = cell(maxk,1);
            A = cell(maxk,1);
            for k = 1:maxk
                fprintf('Simulating %d storms from %s(%d)\n',nStorms,method,k);
                % p = 3;
                A{k} = simulateStorms_v22({mdl.SP, mdl.AP.fit, mdl.F.fit, mdl.R.fit}, nStorms, k, data, method, simulationType, rsdModel, peaks);
%                 storms.(method){k} = simulateStorms_v22({mdl.SP, mdl.AP.fit, mdl.F.fit, mdl.R.fit}, nStorms, p, data, method, simulationType, rsdModel, peaks);
            end
            storms.(method) = A;
        otherwise
            k = 1;
%             maxp = 6;
%             storms.(method) = cell(maxp,1);
            A = cell(maxp,1);

            for p = 1:maxp
                fprintf('Simulating %d storms from %s(%d,%d)\n',nStorms,method,p,k);
                A{p} = simulateStorms_v22({mdl.SP, mdl.AP.fit, mdl.F.fit, mdl.R.fit}, nStorms, p, data, method, simulationType, rsdModel, peaks);
%                 storms.(method){p} = simulateStorms_v22({mdl.SP, mdl.AP.fit, mdl.F.fit, mdl.R.fit}, nStorms, k, data, method, simulationType, rsdModel, peaks);
            end
            storms.(method) = A;
    end
end

end