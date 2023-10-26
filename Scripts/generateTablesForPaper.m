%%
MatPrm = cellfun(@(x)(cell2mat(x.AP.fit{4}.phat')),myMdlFits,'UniformOutput',false);
MatX = cell(size(MatPrm));
MatY = cell(size(MatPrm));
for i = 1:length(MatPrm)
    nPrm = size(MatPrm{i},2);
    nX = floor(nPrm/2);
    MatX{i} = [MatPrm{i}(:,1:nX/2),[1;-inf],MatPrm{i}(:,nX/2+1:nX)];
    MatY{i} = MatPrm{i}(:,nX+1:nPrm);
end

fprintf('\nAlpha matrix:\n');
fprintf('\\begin{pmatrix}\n');
for iPrm = 1:nX+1
    
    x = cellfun(@(x)(x(1,iPrm)),MatX);
    y = cellfun(@(x)(x(1,iPrm)),MatY);
    
    if iPrm == nX/2+1
        fprintf('& %.2f (%.2f,%.2f) \\\\\n',y(1),quantile(y,0.05),quantile(y,0.95));
    else
        fprintf('%.2f (%.2f,%.2f) & %.2f (%.2f,%.2f) \\\\\n',x(1),quantile(x,0.05),quantile(x,0.95),y(1),quantile(y,0.05),quantile(y,0.95));
    end
end
fprintf('\\end{pmatrix}\n');
fprintf('\n\nBeta matrix:\n');
fprintf('\\begin{pmatrix}\n');
for iPrm = 1:nX+1
    
    x = cellfun(@(x)(x(2,iPrm)),MatX);
    y = cellfun(@(x)(x(2,iPrm)),MatY);
    
    if iPrm == nX/2+1
        fprintf('& %.2f (%.2f,%.2f) \\\\\n',y(1),quantile(y,0.05),quantile(y,0.95));
    else
        fprintf('%.2f (%.2f,%.2f) & %.2f (%.2f,%.2f) \\\\\n',x(1),quantile(x,0.05),quantile(x,0.95),y(1),quantile(y,0.05),quantile(y,0.95));
    end
end
fprintf('\\end{pmatrix}\n');

%%
% MMEM
mleX = myMdlFits{1}.(s).fit.typen{2}.X{1}.newMLE;
mleY = myMdlFits{1}.(s).fit.typen{2}.Y{1}.newMLE;

mleX = cellfun(@(y)(cell2mat(cellfun(@(x)(x.newMLE),y.(s).fit.type1.X,'UniformOutput',false)')),myMdlFits,'UniformOutput',false);
mleY = cellfun(@(y)(cell2mat(cellfun(@(x)(x.newMLE),y.(s).fit.type1.Y,'UniformOutput',false)')),myMdlFits,'UniformOutput',false);

fprintf('\nAlpha matrix \n\\begin{tabular}{c c}\n');
for i = 1:5
    currY = cellfun(@(y)(y(1,i)),mleY);
    if i == 1
        fprintf('  & %.2f (%.2f, %.2f) \\\\\n',currY(1),quantile(currY,0.05),quantile(currY,0.95));
    else
        currX = cellfun(@(y)(y(1,i-1)),mleX);
        fprintf(' %.2f (%.2f, %.2f) & %.2f (%.2f, %.2f) \\\\\n',currX(1),quantile(currX,0.05),quantile(currX,0.95),currY(1),quantile(currY,0.05),quantile(currY,0.95));
    end
end
fprintf('\\end{tabular}\n');

fprintf('\n\nBeta Matrix\n\\begin{tabular}{c c}\n');
for i = 1:5
    currY = cellfun(@(y)(y(2,i)),mleY);
    if i == 1
        fprintf(' & %.2f (%.2f, %.2f) \\\\\n',currY(1),quantile(currY,0.05),quantile(currY,0.95));
    else
        currX = cellfun(@(y)(y(2,i-1)),mleX);
        fprintf(' %.2f (%.2f, %.2f) & %.2f (%.2f, %.2f) \\\\\n',currX(1),quantile(currX,0.05),quantile(currX,0.95),currY(1),quantile(currY,0.05),quantile(currY,0.95));
    end
end
fprintf('\\end{tabular}\n');

%% EVAR1
EVAR1_X =  cellfun(@(x)(x.(s).fit.typen{1}.X{1}.newMLE),myMdlFits,'UniformOutput',false);
EVAR1_Y =  cellfun(@(x)(x.(s).fit.typen{1}.Y{1}.newMLE),myMdlFits,'UniformOutput',false);
EVAR1_B = cell2mat(cellfun(@(x,y)([x(end),y(end)]),EVAR1_X,EVAR1_Y,'UniformOutput',false));
EVAR1_Phi1 = cellfun(@(x,y)([x(1),x(2);y(1),y(2)]),EVAR1_X,EVAR1_Y,'UniformOutput',false);
fprintf('\n\nEVAR(1) Matrix\n\\begin{tabular}{c | c c}\n');
fprintf('EVAR(1) \\\\ \\hline\n')
for i = 1:2
    if i == 1
        fprintf('$\\Phi^{(1)}$ &');
    else
        fprintf(' &');
    end
    curr1 = cellfun(@(y)(y(i,1)),EVAR1_Phi1);
    curr2 = cellfun(@(y)(y(i,2)),EVAR1_Phi1);
    fprintf('%.2f (%.2f, %.2f) & %.2f (%.2f, %.2f) \\\\\n',curr1(1),quantile(curr1,0.05),quantile(curr1,0.95),curr2(1),quantile(curr2,0.05),quantile(curr2,0.95));
end
fprintf('\\hline\n')
fprintf('$\\bold{B}$ &');
curr1 = EVAR1_B(:,1);
curr2 = EVAR1_B(:,2);
fprintf('%.2f (%.2f, %.2f) & %.2f (%.2f, %.2f) \\\\\n',curr1(1),quantile(curr1,0.05),quantile(curr1,0.95),curr2(1),quantile(curr2,0.05),quantile(curr2,0.95));
fprintf('\\end{tabular}\n');


%% EVAR2
% EVAR2_X =  myMdlFits{1}.(s).fit.typen{2}.X{1}.newMLE;
% EVAR2_Y =  myMdlFits{1}.(s).fit.typen{2}.Y{1}.newMLE;
EVAR2_X =  cellfun(@(x)(x.(s).fit.typen{2}.X{1}.newMLE),myMdlFits,'UniformOutput',false);
EVAR2_Y =  cellfun(@(x)(x.(s).fit.typen{2}.Y{1}.newMLE),myMdlFits,'UniformOutput',false);

% EVAR2_B = [EVAR2_X(end),EVAR2_Y(end)];
EVAR2_B = cell2mat(cellfun(@(x,y)([x(end),y(end)]),EVAR2_X,EVAR2_Y,'UniformOutput',false));

EVAR2_Phi1 = cellfun(@(x,y)([x(2),x(4);y(2),y(4)]),EVAR2_X,EVAR2_Y,'UniformOutput',false);
EVAR2_Phi2 = cellfun(@(x,y)([x(1),x(3);y(1),y(3)]),EVAR2_X,EVAR2_Y,'UniformOutput',false);

fprintf('\n\nEVAR(2) Matrix\n\\begin{tabular}{c | c c}\n');
fprintf('EVAR(2) \\\\ \\hline\n')
for i = 1:2
    if i == 1
        fprintf('$\\Phi^{(1)}$ &');
    else
        fprintf(' &');
    end
    curr1 = cellfun(@(y)(y(i,1)),EVAR2_Phi1);
    curr2 = cellfun(@(y)(y(i,2)),EVAR2_Phi1);
    fprintf('%.2f (%.2f, %.2f) & %.2f (%.2f, %.2f) \\\\\n',curr1(1),quantile(curr1,0.05),quantile(curr1,0.95),curr2(1),quantile(curr2,0.05),quantile(curr2,0.95));
end
fprintf('\\hline\n')
for i = 1:2
    if i == 1
        fprintf('$\\Phi^{(2)}$ &');
    else
        fprintf(' &');
    end
    curr1 = cellfun(@(y)(y(i,1)),EVAR2_Phi2);
    curr2 = cellfun(@(y)(y(i,2)),EVAR2_Phi2);
    fprintf('%.2f (%.2f, %.2f) & %.2f (%.2f, %.2f) \\\\\n',curr1(1),quantile(curr1,0.05),quantile(curr1,0.95),curr2(1),quantile(curr2,0.05),quantile(curr2,0.95));
end
fprintf('\\hline\n')
fprintf('$\\bold{B}$ &');
curr1 = EVAR2_B(:,1);
curr2 = EVAR2_B(:,2);
fprintf('%.2f (%.2f, %.2f) & %.2f (%.2f, %.2f) \\\\\n',curr1(1),quantile(curr1,0.05),quantile(curr1,0.95),curr2(1),quantile(curr2,0.05),quantile(curr2,0.95));
fprintf('\\end{tabular}\n');

%% EVAR4
evarOrd = 4;
EVAR_X =  cellfun(@(x)(x.(s).fit.typen{evarOrd}.X{1}.newMLE),myMdlFits,'UniformOutput',false);
EVAR_Y =  cellfun(@(x)(x.(s).fit.typen{evarOrd}.Y{1}.newMLE),myMdlFits,'UniformOutput',false);

% EVAR2_B = [EVAR2_X(end),EVAR2_Y(end)];
EVAR_B = cell2mat(cellfun(@(x,y)([x(end),y(end)]),EVAR_X,EVAR_Y,'UniformOutput',false));

EVAR_Phi1 = cellfun(@(x,y)([x(4),x(8);y(4),y(8)]),EVAR_X,EVAR_Y,'UniformOutput',false);
EVAR_Phi2 = cellfun(@(x,y)([x(3),x(7);y(3),y(7)]),EVAR_X,EVAR_Y,'UniformOutput',false);
EVAR_Phi3 = cellfun(@(x,y)([x(2),x(6);y(2),y(6)]),EVAR_X,EVAR_Y,'UniformOutput',false);
EVAR_Phi4 = cellfun(@(x,y)([x(1),x(5);y(1),y(5)]),EVAR_X,EVAR_Y,'UniformOutput',false);

fprintf('\n\nEVAR(%d) Matrix\n\\begin{tabular}{c | c c}\n',evarOrd);
fprintf('EVAR(%d) \\\\ \\hline\n',evarOrd)
for i = 1:2
    if i == 1
        fprintf('$\\Phi^{(1)}$ &');
    else
        fprintf(' &');
    end
    curr1 = cellfun(@(y)(y(i,1)),EVAR_Phi1);
    curr2 = cellfun(@(y)(y(i,2)),EVAR_Phi1);
    fprintf('%.2f (%.2f, %.2f) & %.2f (%.2f, %.2f) \\\\\n',curr1(1),quantile(curr1,0.05),quantile(curr1,0.95),curr2(1),quantile(curr2,0.05),quantile(curr2,0.95));
end
fprintf('\\hline\n')
for i = 1:2
    if i == 1
        fprintf('$\\Phi^{(2)}$ &');
    else
        fprintf(' &');
    end
    curr1 = cellfun(@(y)(y(i,1)),EVAR_Phi2);
    curr2 = cellfun(@(y)(y(i,2)),EVAR_Phi2);
    fprintf('%.2f (%.2f, %.2f) & %.2f (%.2f, %.2f) \\\\\n',curr1(1),quantile(curr1,0.05),quantile(curr1,0.95),curr2(1),quantile(curr2,0.05),quantile(curr2,0.95));
end
fprintf('\\hline\n')
for i = 1:2
    if i == 1
        fprintf('$\\Phi^{(3)}$ &');
    else
        fprintf(' &');
    end
    curr1 = cellfun(@(y)(y(i,1)),EVAR_Phi3);
    curr2 = cellfun(@(y)(y(i,2)),EVAR_Phi3);
    fprintf('%.2f (%.2f, %.2f) & %.2f (%.2f, %.2f) \\\\\n',curr1(1),quantile(curr1,0.05),quantile(curr1,0.95),curr2(1),quantile(curr2,0.05),quantile(curr2,0.95));
end
fprintf('\\hline\n')
for i = 1:2
    if i == 1
        fprintf('$\\Phi^{(4)}$ &');
    else
        fprintf(' &');
    end
    curr1 = cellfun(@(y)(y(i,1)),EVAR_Phi4);
    curr2 = cellfun(@(y)(y(i,2)),EVAR_Phi4);
    fprintf('%.2f (%.2f, %.2f) & %.2f (%.2f, %.2f) \\\\\n',curr1(1),quantile(curr1,0.05),quantile(curr1,0.95),curr2(1),quantile(curr2,0.05),quantile(curr2,0.95));
end
fprintf('\\hline\n')
fprintf('$\\bold{B}$ &');
curr1 = EVAR_B(:,1);
curr2 = EVAR_B(:,2);
fprintf('%.2f (%.2f, %.2f) & %.2f (%.2f, %.2f) \\\\\n',curr1(1),quantile(curr1,0.05),quantile(curr1,0.95),curr2(1),quantile(curr2,0.05),quantile(curr2,0.95));
fprintf('\\end{tabular}\n');

