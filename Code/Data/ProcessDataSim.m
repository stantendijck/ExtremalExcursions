function data = ProcessDataSim(stormsOrg, thr)
% Get storms

% Storm threshold
if nargin < 2
    q = 0.95;
end


methods = fieldnames(stormsOrg);
% nOrder = size(stormsOrg.(methods{1}),1);


for iMethod = 1:length(methods)
    nOrder = length(stormsOrg.(methods{iMethod}));
	data.storms.(methods{iMethod}) = cell(nOrder,1);
	for iOrder = 1:nOrder
		theseStorms = stormsOrg.(methods{iMethod}){iOrder};
		storms = {};
		for iStorm = 1:length(theseStorms)
			Y = theseStorms{iStorm}.originalMargins;
			D = theseStorms{iStorm}.direction;
            n = size(Y,1);
			
			% Initialisation
			% stormIn = data.lap.X(1) > data.thr | data.lap.Y(1) > data.thr;
			stormIn = Y(1,1) > thr;
			if stormIn
				storm.data = [Y(1,1),Y(1,2)];
                storm.direction = [D(1,1),D(1,2)];
				storm.index = 1;
			else
				storm.data = nan(0,2);
                storm.direction = nan(0,2);
				storm.index = [];
			end


			for i = 2:n
				%if stormIn && (data.lap.X(i) > thr || data.lap.Y(i) > thr)
				if stormIn && Y(i,1) > thr
					storm.data(end+1,:) = Y(i,:);
					storm.direction(end+1,:) = D(i,:);
					storm.index(end+1) = i;
				elseif stormIn
					storms{end+1} = storm;
					stormIn = false;
				%elseif data.lap.X(i) > thr || data.lap.Y(i) > thr
				elseif Y(i,1) > thr
					storm.data = Y(i,:);
					storm.direction = D(i,:);
					storm.index = i;
					stormIn = true;
				end
				%if (data.lap.X(i) > thr || data.lap.Y(i) > thr) && i == n
				if Y(i,1) > thr && i == n
					storms{end+1} = storm;
				end
			end
		end
		data.storms.(methods{iMethod}){iOrder} = storms;
	end
end

end