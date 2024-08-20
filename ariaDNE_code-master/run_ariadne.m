
path = '/Users/rileywilde/SAMS/teeth_sep22/Default/HecateVertex/NumSegments_5/NumEigs_2/segment';

addpath('./utils')
addpath(genpath('./utils'))

addpath(path)
addpath(genpath(path))

%[fileNameList, suffix] = getFileNames('./alginedMeshes2');

addpath('./alginedMeshes2')
%ddd = dir('./alignedMeshes2')

[fileNameList,suffix] = getFileNames(['alignedMeshes2'])

Options.distInfo = 'Geodeisic';
Options.cutThresh = 0;
bandwidth = 0.06;




seg = [{'_seg1'},{'_seg2'},{'_seg3'},{'_seg4'},{'_seg5'}];
result = zeros(length(fileNameList),length(seg));

for j = 1:5
    s = seg{j};
    for i = 1:length(fileNameList)
        disp([j,i])
        meshname = [path,'/',fileNameList{i}, s, suffix];
        H = ariaDNE(meshname, bandwidth, Options);
        result(i,j) = H.dne;
    end
end


header = cat(2,{'fname'},seg);
C = cat(2,fileNameList',num2cell(result));

T = cell2table(C,'VariableNames',header);
writetable(T,'ariadne_teeth_segs.csv')

