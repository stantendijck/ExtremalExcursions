function mdl = fitModels_v3(htdata, printOn)

%% Initialisation
if nargin < 2
    printOn = true;
end
adpMCMCFlag = false;

%% Model of type 1
data.explanatory = htdata.X.data(:,1);
data.logExplanatory = log(data.explanatory);

% Fit Y(i) | X(1) > u                                   (i = 1, 2, ...)
mdl = struct();
mdl.type1 = struct();
mdl.type1.Y = fitMdl1(data, htdata.Y.data, adpMCMCFlag, {printOn, 'Y', 0});

% Fit X(i+1) | X(1) > u                                 (i = 1, 2, ...)
mdl.type1.X = fitMdl1(data, htdata.X.data(:,2:end), adpMCMCFlag, {printOn, 'X', 1});


%% Model of type n
maxType = size(htdata.Y.data,2)-1;
mdl.typen = cell(maxType-1,1);
mdl.varn = cell(maxType-1,1);
for type = 2:maxType
    data.explanatory = [htdata.X.data(:,1:type-1), htdata.Y.data(:,1:type-1)];
    data.logExplanatory = log(data.explanatory(:,1));
    
    %ind = type:size(htdata.X.data,2);
    ind = type;
 
    % Fit X(i+1) | (X(1) > u, Y(1),..., X(i), Y(i))         (i = 1, 2, ...)
    mdl.varn{type-1} = struct();
    mdl.varn{type-1}.X = fitVAR_v2(data, htdata.X.data(:,ind), mdl.type1, 'X', type, adpMCMCFlag, {printOn, 'X', type-1});

    % Fit Y(i+1) | (X(1) > u, Y(1),..., X(i), Y(i))         (i = 1, 2, ...)
    mdl.varn{type-1}.Y = fitVAR_v2(data, htdata.Y.data(:,ind), mdl.type1, 'Y', type, adpMCMCFlag, {printOn, 'Y', type-1});

    % Fit X(i+1) | (X(1) > u, Y(1),..., X(i), Y(i))         (i = 1, 2, ...)
    mdl.typen{type-1} = struct();
    mdl.typen{type-1}.X = fitMdln_v2(data, htdata.X.data(:,ind), mdl.type1, mdl.varn{type-1}.X, 'X', type, adpMCMCFlag, {printOn, 'X', type-1});

    % Fit Y(i+1) | (X(1) > u, Y(1),..., X(i), Y(i))         (i = 1, 2, ...)
    mdl.typen{type-1}.Y = fitMdln_v2(data, htdata.Y.data(:,ind), mdl.type1, mdl.varn{type-1}.Y, 'Y', type, adpMCMCFlag, {printOn, 'Y', type-1});
    
end


end

