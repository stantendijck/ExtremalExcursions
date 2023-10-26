function keefProp = generateKeefProp(respData, explData)

keefProp = struct();
keefProp.plus.data = respData - explData(:,1);
keefProp.min.data = respData + explData(:,1);
keefProp.plus.max = max(keefProp.plus.data);
keefProp.plus.min = min(keefProp.plus.data);
keefProp.min.max = max(keefProp.min.data);
keefProp.min.min = min(keefProp.min.data);
keefProp.X.max = max(explData(:,1));

end