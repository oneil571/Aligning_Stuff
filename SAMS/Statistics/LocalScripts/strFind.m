function ind = strFind(str,dict)
for i = 1:length(dict)
    if strcmp(str,dict{i})
        ind = i;
        break;
    elseif i == length(dict)
        ind = -1;
        disp([str ' not found in dictionary']);
    end
end
end