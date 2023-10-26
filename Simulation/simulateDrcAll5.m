function stormsOrg = simulateDrcAll5(storms, mdldrc, data, MMX, MMY, I)

stormsOrg = struct();
fnames = fieldnames(storms);
upcrossingInds = cellfun(@(x)(x.index(1)),data.storms);
for iName = 1:length(fnames)
    fprintf('%s\n',fnames{iName});
    maxOrder = length(storms.(fnames{iName}));
    stormsOrg.(fnames{iName}) = cell(maxOrder,1);
    for iOrder = 1:maxOrder
        fprintf('%d/%d\n',iOrder,maxOrder);
        nStorms = length(storms.(fnames{iName}){iOrder});
        stormsOrg.(fnames{iName}){iOrder} = cell(nStorms,1);
%         nStorms = 100;
        currStorms = storms.(fnames{iName}){iOrder};
        currStormsNew = cell(nStorms,1);
        parfor iStorm = 1:nStorms
            if mod(iStorm,100) == 0
                fprintf('+');
            end
            storm = currStorms{iStorm};
            currStormsNew{iStorm} = struct();
            currStormsNew{iStorm}.standardMargins = storm;
            drc = simulateDrc5(storm, mdldrc, data, MMX, MMY, upcrossingInds, I(iStorm));
            currStormsNew{iStorm}.originalMargins = drc.data;
            currStormsNew{iStorm}.direction = drc.drcdata;
        end
        stormsOrg.(fnames{iName}){iOrder} = currStormsNew;
        fprintf('\n');
    end
end

end