
function PrintPseudolandmarks(workingPath,numGPLmks,numFPSLmks)

load([workingPath 'newMeshList.mat']);
load([workingPath 'FinalDists.mat']);

[~,frechInd] = min(sum(dists.^2));

frechMean = newMeshList{frechMean};

%% Get random sample based on conf max
numPseudolmks = numGPLmks+numFPSLmks;

FPSInds = unique(frechMean.GetGPLmk(numGPLmks));
FPSInds = frechMean.GeodesicFarthestPointSampling(numGPLmks+numFPSLmks,FPSInds);
% Build CSV
touch([workingPath 'PseudolandmarkFiles']);
pseudoPath = [workingPath 'PseudolandmarkFiles/'];
touch([pseudoPath 'CSV/']); csvPath = [pseudoPath 'CSV/']);
touch([pseudoPath 'Morphologika/']); morphoPath = [pseudoPath 'Morphologika/']);
load([workingPath 'newMeshListNames.mat']);

fid = fopen([csvPath num2str(numGPLmks) '_GP_' num2str(numFPSLmks) '_FPS.csv'],'w');

curLine = ['SpecimenName,'];
for i = 1:numPseudolmks
    curLine = [curLine 'X ' num2str(i) ',Y ' num2str(i) ',Z ' num2str(i)];
    if i == numPseudolmks
        curLine = [curLine '\n'];
    else
        curLine = [curLine ','];
    end
end
fprintf(fid,curLine);

for n = 1:length(Names)
    curLine = [Names{n} ','];
    for i = 1:numPseudolmks
        curLine = [curLine num2str(newMeshList{n}.V(1,FPSInds(i)))];
        curLine = [curLine ',' num2str(newMeshList{n}.V(2,FPSInds(i)))];
        curLine = [curLine ',' num2str(newMeshList{n}.V(3,FPSInds(i)))];
        if i == numPseudolmks
            curLine = [curLine '\n'];
        else
            curLine = [curLine ','];
        end
    end
    fprintf(fid,curLine);
end
fclose(fid);   


% And now morphologika
fid = fopen([morphoPath num2str(numGPLmks) '_GP_' num2str(numFPSLmks) '_FPS.csv'],'w');
fprintf(fid,'[Individuals]\n');
fprintf(fid,[num2str(length(Names)) '\n']);
fprintf(fid,'[landmarks]\n');
fprintf(fid,[num2str(numGPLmks+numFPSLmks) '\n']);
fprintf(fid,'[dimensions]\n3\n[names]\n');
for n = 1:length(Names)
    fprintf(fid,[Names{n} '\n']);
end
fprintf(fid,'\n[rawpoints]\n\n');
for n = 1:length(Names)
    fprintf(fid,['''' Names{n} '\n\n']);
    for i = 1:numPseudolmks
        curStr = [num2str(newMeshList{n}.V(1,FPSInds(i))) ' ' ...
            num2str(newMeshList{n}.V(2,FPSInds(i))) ' ' ...
            num2str(newMeshList{n}.V(3,FPSInds(i))) '\n'];
        fprintf(fid,curStr);
        if i == numPseudolmks && n < length(Names)
            fprintf(fid,'\n');
        end
    end
end
fclose(fid);