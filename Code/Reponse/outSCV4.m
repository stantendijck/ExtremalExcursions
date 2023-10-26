function impact = outSCV4(testData, stormsOrg, trainData, setZero)

%% Storm trajectory

% xL = linspace(0,360,1000);
% yL = arrayfun(@(x)(transform2Org(Laplace_iCDF(0.95), x, MM_HS)),xL);

% Let's say storm has impact if Hs > 10 and additional impact if WS > 20
% impact:
% Hs < HsThr:    Impact = Inline_WS^2
% Hs > HsThr:    Impact = A * AngHSdrc * (Hs - HsThr) * Hs^2 + Inline_WS^2

%%
% MdlPrm = {[0.4, 0.2], [5, 7]};
% MdlPrm = {[0.3, 0.15], [6, 8]};
MdlPrm = {[0.2, 0.3], [6, 7]};
funs = {@max, @sum};

% impact = struct;
if nargin < 4
    setZero = -1;
end

impact = cell(length(MdlPrm),2);
for iPrm = 1:length(MdlPrm)
    for iFun = 1:2
        % Pre Allocation
        impact{iPrm,iFun} = struct;
        %         impact{iPrm,iFun}.test = ;
        %         impact{iPrm,iFun}.train = ;
        %         impact{iPrm,iFun}.models = ;
        %         impact{iPrm,iFun}.HM = ;
        
        % Response parameters
        A = MdlPrm{iPrm}(1);
        HsThr = MdlPrm{iPrm}(2);
        fun = funs{iFun};
        
        
        %%
        methods = fieldnames(stormsOrg);
        maxModelOrder = 6;
        impact{iPrm,iFun}.models = cell(length(methods) - 1, maxModelOrder);
        impact{iPrm,iFun}.HM = cell(1,1);
        figcnt = 1;
        for iMethod = 1:length(methods)
            if iMethod <= 3
                maxModelOrder = length(stormsOrg.(methods{iMethod}));
                for iOrder = 1:maxModelOrder
                    nStorms = length(stormsOrg.(methods{iMethod}){iOrder});
                    impact{iPrm,iFun}.models{iMethod,iOrder} = nan(nStorms,1);
                    for iStorm = 1:nStorms
                        currStorm = stormsOrg.(methods{iMethod}){iOrder}{iStorm}.originalMargins;
                        %     currStorm = [data.X(data.storms{iStorm}.index),data.Y(data.storms{iStorm}.index)];
                        I = currStorm(:,1) > HsThr & currStorm(:,1) < max(currStorm(:,1));
                        [~,maxInd] = max(currStorm(:,1));
                        I(max(1,maxInd-setZero):min(end,maxInd+setZero)) = 0;
                        notI = ~I;
                        notI(max(1,maxInd-setZero):min(end,maxInd+setZero)) = 0;
                        
                        
                        currStormDrc = stormsOrg.(methods{iMethod}){iOrder}{iStorm}.direction;
                        R = zeros(size(currStorm,1),1);
                        c = cos(abs(currStormDrc(:,2)-currStormDrc(:,1)) * 2 * pi / 360);
                        inline_windspeed = c .* currStorm(:,2);

                        if any(I)
                            AngDrc = 1./cos((abs(mod(stormsOrg.(methods{iMethod}){iOrder}{iStorm}.direction(I,1)+45,90)-45))/360 * 2 * pi);
                            % c = cos(abs(stormsOrg.(methods{iMethod}){iOrder}{iStorm}.direction(I,2)-stormsOrg.(methods{iMethod}){iOrder}{iStorm}.direction(I,1)) * 2 * pi / 360);
                            %inline_windspeed = c .* currStorm(I,2);
                            % [impact(iStorm,iMethod,iOrder),ind] = max(A * AngDrc .* (currStorm(I,1) - HsThr) .* currStorm(I,1).^2 + inline_windspeed.^2);
                            %currImpactI = fun(AngDrc .* (currStorm(I,1) - HsThr) .* currStorm(I,1).^2 + A * inline_windspeed.^2);
                            R(I) = AngDrc .* (currStorm(I,1) - HsThr) .* currStorm(I,1).^2 + A * inline_windspeed(I).^2;
                        end
                        if any(notI)
                            % currImpactnotI = fun(currStorm(~I,2).^2);
                            R(notI) = A * inline_windspeed(notI).^2;
                        end
                        %currImpact = fun([currImpactI;currImpactnotI]);
                        currImpact = fun(R);
                        impact{iPrm,iFun}.models{iMethod,iOrder}(iStorm) = currImpact;
%                         else
%                             % impact(iStorm,iMethod,iOrder) = 0; %sum(C * currStorm(:,2));
%                             currImpactnotI = fun(A * currStorm(~I,2).^2);
% 
%                             impact{iPrm,iFun}.models{iMethod,iOrder}(iStorm) = currImpactnotI;
%                         end
                    end
                    figcnt = figcnt + 1;
                end
            else
                nStorms = length(stormsOrg.HM{1});
                impact{iPrm,iFun}.HM{1} = nan(nStorms,1);
                for iStorm = 1:nStorms
                    currStorm = stormsOrg.HM{1}{iStorm}.originalMargins;
                    I = currStorm(:,1) > HsThr;
                    [~,maxInd] = max(currStorm(:,1));
                    I(max(1,maxInd-setZero):min(end,maxInd+setZero)) = 0;
                    notI = ~I;
                    notI(max(1,maxInd-setZero):min(end,maxInd+setZero)) = 0;
                    
                    currStormDrc = stormsOrg.HM{1}{iStorm}.direction;
                    R = zeros(size(currStorm,1),1);
                    c = cos(abs(currStormDrc(:,2)-currStormDrc(:,1)) * 2 * pi / 360);
                    inline_windspeed = c .* currStorm(:,2);
                    
                    if any(I)
                        AngDrc = 1./cos((abs(mod(stormsOrg.HM{1}{iStorm}.direction(I,1)+45,90)-45))/360 * 2 * pi);
                        %c = cos(abs(stormsOrg.HM{1}{iStorm}.direction(I,2)-stormsOrg.HM{1}{iStorm}.direction(I,1)) * 2 * pi / 360);
                        %inline_windspeed = c .* currStorm(I,2);
                        %currImpactI = fun(AngDrc .* (currStorm(I,1) - HsThr) .* currStorm(I,1).^2 + A * inline_windspeed.^2);
                        R(I) = AngDrc .* (currStorm(I,1) - HsThr) .* currStorm(I,1).^2 + A * inline_windspeed(I).^2;
                    end
                    if any(notI)
                        %currImpactnotI = fun(currStorm(~I,2).^2);
                        R(notI) = A * inline_windspeed(notI).^2;
                    end
                    %currImpact = fun([currImpactI;currImpactnotI]);
                    currImpact = fun(R);
                        
                    impact{iPrm,iFun}.HM{1}(iStorm) = currImpact;
                    %else
                    %    currImpactnotI = fun(A * currStorm(~I,2).^2);
                    %
                    %    impact{iPrm,iFun}.HM{1}(iStorm) = currImpactnotI; %sum(C * currStorm(:,2));
                    %end
                end
            end
        end
        
        nStorms = length(testData.storms);
        impact{iPrm,iFun}.test = cell(1,1);
        impact{iPrm,iFun}.test{1} = nan(nStorms,1);
        for iStorm = 1:nStorms
            currStorm = [testData.X(testData.storms{iStorm}.index),testData.Y(testData.storms{iStorm}.index)];
            I = currStorm(:,1) > HsThr;
            [~,maxInd] = max(currStorm(:,1));
            I(max(1,maxInd-setZero):min(end,maxInd+setZero)) = 0;
            notI = ~I;
            notI(max(1,maxInd-setZero):min(end,maxInd+setZero)) = 0;
            
            currStormDrc = [testData.Xdrc(testData.storms{iStorm}.index),testData.Ydrc(testData.storms{iStorm}.index)];
            R = zeros(size(currStorm,1),1);
            c = cos(abs(currStormDrc(:,2)-currStormDrc(:,1)) * 2 * pi / 360);
            inline_windspeed = c .* currStorm(:,2);

            if any(I)
                AngDrc = 1./cos((abs(mod(currStormDrc(I,1)+45,90)-45))/360 * 2 * pi);
                %c = cos(abs(currStormDrc(I,2)-currStormDrc(I,1)) * 2 * pi / 360);
                %inline_windspeed = c .* currStorm(I,2);
                %currImpactI = fun(AngDrc .* (currStorm(I,1) - HsThr) .* currStorm(I,1).^2 + A * inline_windspeed.^2);
                R(I) = AngDrc .* (currStorm(I,1) - HsThr) .* currStorm(I,1).^2 + A * inline_windspeed(I).^2;
            end
            if any(notI)
                %currImpactnotI = fun(currStorm(~I,2).^2);
                R(notI) = A * inline_windspeed(notI).^2;
            end
            %currImpact = fun([currImpactI;currImpactnotI]);
            currImpact = fun(R);
                
            impact{iPrm,iFun}.test{1}(iStorm) = currImpact;
            %else
            %    currImpactnotI = fun(A * currStorm(~I,2).^2);
            %
            %    impact{iPrm,iFun}.test{1}(iStorm) = currImpactnotI; %sum(C*currStorm(:,2));
            %end
        end
        
        if nargin >= 3
            nStorms = length(trainData.storms);
            impact{iPrm,iFun}.train = cell(1,1);
            impact{iPrm,iFun}.train{1} = nan(nStorms,1);
            for iStorm = 1:nStorms
                currStorm = [trainData.X(trainData.storms{iStorm}.index),trainData.Y(trainData.storms{iStorm}.index)];
                I = currStorm(:,1) > HsThr;
                [~,maxInd] = max(currStorm(:,1));
                I(max(1,maxInd-setZero):min(end,maxInd+setZero)) = 0;
                notI = ~I;
                notI(max(1,maxInd-setZero):min(end,maxInd+setZero)) = 0;
                
                currStormDrc = [trainData.Xdrc(trainData.storms{iStorm}.index),trainData.Ydrc(trainData.storms{iStorm}.index)];
                R = zeros(size(currStorm,1),1);
                c = cos(abs(currStormDrc(:,2)-currStormDrc(:,1)) * 2 * pi / 360);
                inline_windspeed = c .* currStorm(:,2);

                if any(I)
                    AngDrc = 1./cos((abs(mod(currStormDrc(I,1)+45,90)-45))/360 * 2 * pi);
                    %c = cos(abs(currStormDrc(I,2)-currStormDrc(I,1)) * 2 * pi / 360);
                    %inline_windspeed = c .* currStorm(I,2);
                    %currImpactI = fun(AngDrc .* (currStorm(I,1) - HsThr) .* currStorm(I,1).^2 + A * inline_windspeed.^2);
                    R(I) = AngDrc .* (currStorm(I,1) - HsThr) .* currStorm(I,1).^2 + A * inline_windspeed(I).^2;
                end
                if any(notI)
                    %currImpactnotI = fun(currStorm(~I,2).^2);
                    R(notI) = A * inline_windspeed(notI).^2;
                end
                %currImpact = fun([currImpactI;currImpactnotI]);
                currImpact = fun(R);
                    
                impact{iPrm,iFun}.train{1}(iStorm) = currImpact;
                %else
                %    currImpactnotI = fun(A * currStorm(~I,2).^2);
                %
                %    impact{iPrm,iFun}.train{1}(iStorm) = currImpactnotI; %sum(C*currStorm(:,2));
                %end
            end
        end
    end
end



end



