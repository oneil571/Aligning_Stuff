%% Global analysis using Weighted Shape Means

%% Load meshes and normalize
meshCollection = cell(1,length(SpecimenTypes));
for i = 1:length(SpecimenTypes)
    load([projectDir SpecimenTypes{i} '/newMeshList.mat']);
    for j = 1:length(newMeshList)
        newMeshList{j}.V = newMeshList{j}.V - ...
            repmat(mean(newMeshList{j}.V')',1,newMeshList{j}.nV);
        newMeshList{j}.V = newMeshList{j}.V/norm(newMeshList{j}.V,'fro');
    end
    meshCollection{i} = newMeshList;
end
%% Compute means for each possible type of key label combination
shapeMeans = cell(length(keys),length(SpecimenTypes));
for k = 1:length(keys)
    currentKey = keys{k};
    totalUnLabels = {};
    for i = 1:length(SpecimenTypes)
        curMeshes = meshCollection{i};
        load([projectDir SpecimenTypes{i} '/Groups.mat']);
        totalLabels{k,i} = lower(Groups(currentKey));
        curLabels = totalLabels{k,i};
        totalUnLabels = unique([totalUnLabels;unique(totalLabels{k,i})]);
        curUnLabels = unique(totalLabels{k,i});
        curAveMeanVerts = zeros(size(newMeshList{1}.V));
        for j = 1:length(curUnLabels)
            inds = find(strcmp(curUnLabels{j},curLabels));
            curMeanVerts = zeros(size(newMeshList{1}.V));
            for p = 1:length(inds)
                curMeanVerts = curMeanVerts + curMeshes{inds(p)}.V;
            end
            curMeanVerts = curMeanVerts/length(inds);
            curMeanVerts = curMeanVerts - ...
                repmat(mean(curMeanVerts')',1,newMeshList{1}.nV);
            curMeanVerts = curMeanVerts/norm(curMeanVerts,'fro');
            curAveMeanVerts = curAveMeanVerts + curMeanVerts;
        end
        curAveMeanVerts = curAveMeanVerts/length(curUnLabels);
        curAveMeanVerts = curAveMeanVerts - ...
            repmat(mean(curAveMeanVerts')',1,newMeshList{1}.nV);
        curAveMeanVerts = curAveMeanVerts/norm(curAveMeanVerts,'fro');
        shapeMeans{k,i} = curAveMeanVerts;
    end
end



dists = cell(length(keys),length(SpecimenTypes));
for i = 1:length(SpecimenTypes)
    for k = 1:length(keys)
        curDists = zeros(length(meshCollection{i}),1);
        for j = 1:length(meshCollection{i})
            curDists(j) = norm(meshCollection{i}{j}.V-shapeMeans{k,i},'fro');
        end
        dists{k,i} = curDists;
    end
end

%% Levene tests for each pair of MetaGroups: prints nothing if only one group

writePath = [interStatPath 'Total/WeightedMeanShape/'];
quadPath = [writePath 'Quad/']; absPath = [writePath 'Abs/'];
touch(quadPath); touch(absPath);
quadid = fopen([writePath 'Quad/LeveneVals.tsv'],'w');
absid = fopen([writePath 'Abs/LeveneVals.tsv'],'w');
quadDict = {'Mean_Grouping\tSpecimen1\tSpecimen2\tP_value\n'}; 
absDict = {'Mean_Grouping\tSpecimen1\tSpecimen2\tP_value\n'}; quadP = []; absP = [];

for k = 1:length(keys)
    for i = 1:length(SpecimenTypes)-1
        for j = i+1:length(SpecimenTypes)
            if isempty(dists{k,i}) || isempty(dists{k,j})
                continue;
            end
            touch([writePath 'Quad/' keys{k}]); 
            touch([writePath 'Abs/' keys{k}]); 
            if length(dists{k,i}) > length(dists{k,j})
                dists_i = dists{k,i};
                dists_j = [dists{k,j};NaN*ones(length(dists{k,i})...
                    -length(dists{k,j}),1)];
                Specimen1 = SpecimenTypes{i};
                Specimen2 = SpecimenTypes{j};
            else
                dists_j = dists{j,m};
                dists_i = [dists{k,i};NaN*ones(length(dists{k,j})...
                    -length(dists{k,i}),1)];
                Specimen1 = SpecimenTypes{j};
                Specimen2 = SpecimenTypes{i};
            end
            P = vartestn([dists_i dists_j],'TestType','LeveneQuadratic',...
            'Display','off');
            save([quadPath keys{k} '/' ...
                SpecimenTypes{i} '_' SpecimenTypes{j} '.mat'],'P');
            quadDict = [quadDict [keys{k} '\t' SpecimenTypes{i} '\t'...
                SpecimenTypes{j} '\t'  num2str(P) '\n']];
            quadP = [quadP P];
            P = vartestn([dists_i dists_j],'TestType','LeveneAbsolute',...
            'Display','off');
            absDict = [absDict [keys{k} '\t' SpecimenTypes{i} '\t'...
                SpecimenTypes{j} '\t'  num2str(P) '\n']];
            absP = [absP P];
            save([absPath keys{k} '/' ...
                SpecimenTypes{i} '_' SpecimenTypes{j} '.mat'],'P');
        end
    end
end

[~,quadInd] = sort(quadP); [~,absInd] = sort(absP);
fprintf(quadid,quadDict{1}); fprintf(absid,absDict{1});
for q = quadInd
    fprintf(quadid,quadDict{q+1});
end
for q = absInd
    fprintf(absid,absDict{q+1});
end
fclose(quadid); fclose(absid);


