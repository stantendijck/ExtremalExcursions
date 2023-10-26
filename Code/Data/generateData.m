function data = generateData(n, modelNumber, stormVars)

switch modelNumber
    case 1
        % Rsd
        [~, x_gumbel, ~] = simulate_alg21(n, 2, 0.5);
        X = nan(n,1); Y = nan(n,1);
        X(1:2) = Exp_iCDF(rand(2,1));
        Y(1:2) = Exp_iCDF(rand(2,1));
        phi = [0.65, 0.3];
        phiTilde = 0.75;
        xi = 0.2;
        for i = 3:n
            X(i) = 0.3 * x_gumbel(i,1) + phi(1) * X(i-1) + phi(2) * X(i-2);
            Y(i) = 0.3 *x_gumbel(i,2) + phiTilde(1) * Y(i-1) + xi(1) * X(i-1);
        end
        
        % Transform data to standard margins
        data.X = X;
        data.Y = Y;
        marg_quant.X.upper = 0.8;
        marg_quant.X.lower = 0.2;
        marg_quant.Y.upper = 0.8;
        marg_quant.Y.lower = 0.2;
        data.unif = PIT(data, marg_quant);
        data.exp.X = Exp_iCDF(data.unif.X);
        data.exp.Y = Exp_iCDF(data.unif.Y);
        
        % Plot
        figure(1);clf;
        plot( data.exp.X, 'k-');
        hold on;
        plot( data.exp.Y, 'r-');
        xlim([0 100]);
        
        % Storm quantile
        q = 0.8;
        thr = Exp_iCDF(q);
        
        % Initialisation
        stormIn = false;
        for stormVar = stormVars
            stormIn = stormIn | data.exp.(stormVar{1})(1);
        end
        %stormIn = data.exp.X(1) > thr | data.exp.Y(1) > thr;
        if stormIn
            storm.data = [data.exp.X(1),data.exp.Y(1)];
            storm.index = 1;
        else
            storm.data = nan(0,2);
            storm.index = [];
        end
        
        % Loop over the data
        storms = {};
        for i = 2:n
            excFlag = 0;
            for stormVar = stormVars
                excFlag = excFlag | data.exp.(stormVar{1})(i) > thr;
            end
            if stormIn && excFlag
                storm.data(end+1,:) = [data.exp.X(i),data.exp.Y(i)];
                storm.index(end+1) = i;
            elseif stormIn
                storms{end+1} = storm;
                stormIn = false;
            elseif excFlag
                storm.data = [data.exp.X(i),data.exp.Y(i)];
                storm.index = i;
                stormIn = true;
            end
            if excFlag && i == n
                storms{end+1} = storm;
            end
        end
        
        figure(2); clf;
        iCnt = 1;
        iStorm = 0;
        while iCnt <= 16
            iStorm = iStorm + 1;
            if length(storms{iStorm}.index)<10
                continue
            end
            subplot(4,4,iCnt);
            inds = storms{iStorm}.index;
            if storms{iStorm}.index(1) > 1
                inds = [inds(1)-1,inds];
            end
            if storms{iStorm}.index(end) < n
                inds = [inds,inds(end)+1];
            end
            plot(data.exp.X(inds));
            hold on;
            plot(data.exp.Y(inds));
            plot(xlim,thr*ones(2,1),'k--');
        	iCnt = iCnt + 1;
        end
        
        data.storms = storms;
        data.thr = thr;
    otherwise
        data = nan;
end

end