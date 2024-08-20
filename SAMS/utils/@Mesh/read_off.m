function [V,F,Fs] = read_off(filename)

fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

str = fgets(fid);   % -1 if eof
% if ~strcmp(str(end-3:end), 'OFF')
if ~findstr(str,'OFF')
    error('The file is not a valid OFF one.');    
end

str = fgets(fid);
[a,str] = strtok(str); Nv= str2num(a);
[a,str] = strtok(str); Nf= str2num(a);



[A,cnt] = fscanf(fid,'%f %f %f', 3*Nv);

if cnt~=3*Nv
    warning('Problem in reading vertices.');
end
V = reshape(A, 3, round(cnt/3));

% read Face 1  1088 480 1022

%% Temporary fix because files are different
testLine = fgets(fid);
testFormat = strsplit(testLine);
if length(testFormat) == 2
    testLine = fgets(fid);
    testFormat = strsplit(testLine);
end

if length(testFormat) == 4 %3 entries
elseif length(testFormat) == 5
else
    error('Problem with mesh');
end
prev = [];
searchStr = '';
len = length(testFormat)-1;
for i = 1:length(testFormat)-1
    prev = [prev;str2num(testFormat{i})];
    searchStr = [searchStr '%d'];
    if i+1 < length(testFormat)
        searchStr = [searchStr ' '];
    else
        searchStr = [searchStr '\n'];
    end
end
[A,cnt] = fscanf(fid,searchStr,len*Nf-len);
A = [prev;A];
if cnt~=(len*Nf-len)
    warning('Problem in reading faces.');
end
F = reshape(A,len,round((cnt+len)/len));
if len == 4
    F = F(2:4,:);
end
F = F-min(min(F))+1;
fclose(fid);


end

