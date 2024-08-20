%% Global analysis using Weighted Shape Means

%% Load MDS information
MDSCollection = cell(1,length(SpecimenTypes));
MeanCoord = cell(1,length(SpecimenTypes));
for i = 1:length(SpecimenTypes)
    load([projectDir SpecimenTypes{i} '/MDSEmbedding.mat']);
    MDSCollection{i} = Y;
end
%% Compute means for each possible type of key label combination
MDSMeans = cell(length(keys),length(SpecimenTypes));
for k = 1:length(keys)
    currentKey = keys{k};
    totalUnLabels = {};
    for i = 1:length(SpecimenTypes)
        curMDS = MDSCollection{i};
        load([projectDir SpecimenTypes{i} '/Groups.mat']);
        totalLabels{k,i} = lower(Groups(currentKey));
        curLabels = totalLabels{k,i};
        totalUnLabels = unique([totalUnLabels;unique(totalLabels{k,i})]);
        curUnLabels = unique(totalLabels{k,i});
        curAveMeanMDS = zeros(1,size(MDSCollection{1},2));
        for j = 1:length(curUnLabels)
            inds = find(strcmp(curUnLabels{j},curLabels));
            curMeanMDS = zeros(1,size(MDSCollection{1},2));
            for p = 1:length(inds)
                curMeanMDS = curMeanMDS + curMDS(inds(p),:);
            end
            curMeanMDS = curMeanMDS/length(inds);
            curAveMeanMDS = curAveMeanMDS + curMeanMDS;
        end
        curAveMeanMDS = curAveMeanMDS/length(curUnLabels);
        MDSMeans{k,i} = curAveMeanMDS;
    end
end



dists = cell(length(keys),length(SpecimenTypes));
for i = 1:length(SpecimenTypes)
    curMDS = MDSCollection{i};
    for k = 1:length(keys)
        curDists = zeros(size(MDSCollection{i},1),1);
        for j = 1:size(MDSCollection{i},1)
            curDists(j) = norm(curMDS(j,:)-MDSMeans{k,i});
        end
        dists{k,i} = curDists;
    end
end

%% Levene tests for each pair of MetaGroups: prints nothing if only one group

writePath = [interStatPath 'Total/WeightedMeanMDS/'];
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


