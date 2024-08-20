%% Global analysis using Sample Mean

%% Load meshes and normalize, then compute distances
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


dists = cell(1,length(SpecimenTypes));
for i = 1:length(SpecimenTypes)
    curMeanVerts = zeros(3,meshCollection{i}{1}.nV);
    for j = 1:length(meshCollection{i})
        curMeanVerts = curMeanVerts + meshCollection{i}{j}.V;
    end
    curMeanVerts = curMeanVerts/length(meshCollection{i});
    %Project back to Procrustes space
    curMeanVerts = curMeanVerts - repmat(mean(curMeanVerts')',1,size(curMeanVerts,2));
    curMeanVerts = curMeanVerts/norm(curMeanVerts,'fro');
    curDists = zeros(length(meshCollection{i}),1);
    for j = 1:length(meshCollection{i})
        curDists(j) = norm(meshCollection{i}{j}.V-curMeanVerts,'fro');
    end
    dists{i} = curDists;
end

%% Levene tests for each pair of MetaGroups: prints nothing if only one group
touch([interStatPath 'Total/SampleMean/Quad/']);
touch([interStatPath 'Total/SampleMean/Abs/']);
quadid = fopen([interStatPath 'Total/SampleMean/Quad/LeveneTests.csv'],'w');
fprintf(quadid,'Specimen1\tSpecimen2\tP_value\n');
absid = fopen([interStatPath 'Total/SampleMean/Abs/LeveneTests.csv'],'w');
fprintf(absid,'Specimen1\tSpecimen2\tP_value\n');
quadDict = {}; absDict = {}; 
quadP = []; absP = [];
for i = 1:length(SpecimenTypes)-1
    for j = (i+1):length(SpecimenTypes)
        if length(dists{i}) > length(dists{j})
            dists_i = dists{i};
            dists_j = [dists{j};NaN*ones(length(dists{i})...
                -length(dists{j}),1)];
            Specimen1 = SpecimenTypes{i};
            Specimen2 = SpecimenTypes{j};
        else
            dists_j = dists{j};
            dists_i = [dists{i};NaN*ones(length(dists{j})...
                -length(dists{i}),1)];
            Specimen1 = SpecimenTypes{j};
            Specimen2 = SpecimenTypes{i};
        end
        P = vartestn([dists_i dists_j],'TestType','LeveneQuadratic',...
            'Display','off');
        quadDict = [quadDict [SpecimenTypes{i} '\t' SpecimenTypes{j} '\t'...
            num2str(P) '\n']];
        quadP = [quadP P];
        touch([interStatPath 'Total/SampleMean/Quad/']);
        save([interStatPath 'Total/SampleMean/Quad/'...
            Specimen1 '_' Specimen2 '.mat'],'P');
        P = vartestn([dists_i dists_j],'TestType','LeveneAbsolute',...
            'Display','off');
        absDict = [absDict [SpecimenTypes{i} '\t' SpecimenTypes{j} '\t'...
            num2str(P) '\n']];
        absP = [absP P];
        touch([interStatPath 'Total/SampleMean/Abs/']);
        save([interStatPath' 'Total/SampleMean/Abs/'...
            Specimen1 '_' Specimen2 '.mat'],'P');
    end
end
[~,quadInd] = sort(quadP); [~,absInd] = sort(absP);
for q = quadInd
    fprintf(quadid,quadDict{q});
end
for q = absInd
    fprintf(absid,absDict{q});
end
fclose(quadid); fclose(absid);
%%
dists = cell(length(keys),length(SpecimenTypes));
totalLabels = cell(length(keys),length(SpecimenTypes));


%% Intergroup analysis for each subgroup
% First step: Obtain different group labels and the lists of unique
% labels
for k = 1:length(keys)
    currentKey = keys{k};
    totalUnLabels = {};
    for i = 1:length(SpecimenTypes)
        load([projectDir SpecimenTypes{i} '/Groups.mat']);
        totalLabels{k,i} = lower(Groups(currentKey));
        totalUnLabels = unique([totalUnLabels;unique(totalLabels{k,i})]);
    end
    %Find dists for all possible combinations of specimens and labels
    dists = cell(length(SpecimenTypes),length(totalUnLabels));
    for i = 1:length(SpecimenTypes)
        for j = 1:length(totalUnLabels)
            inds = find(strcmp(totalUnLabels{j},totalLabels{k,i}));
            if isempty(inds)
                disp(['No ' totalUnLabels{j} ' found in collection of ' ...
                    SpecimenTypes{i}]);
                continue;
            end
            curMeanVerts = zeros(3,meshCollection{i}{1}.nV);
            for m = inds
                curMesh =meshCollection{i}{m};
                curMeanVerts = curMeanVerts + curMesh.V;
            end
            curMeanVerts = curMeanVerts/length(inds);
            %Project back to Procrustes space
            curMeanVerts = curMeanVerts - repmat(mean(curMeanVerts')',1,size(curMeanVerts,2));
            curMeanVerts = curMeanVerts/norm(curMeanVerts,'fro');
            curDists = zeros(length(inds),1);
            for m = 1:length(inds)
                curDists(m) = norm(meshCollection{i}{inds(m)}.V-curMeanVerts,'fro');
            end

            dists{i,j} = curDists;
        end
    end
    %Setup file paths
    writePath = [interStatPath keys{k} '/SampleMean/'];
    touch(writePath);
    touch([writePath 'Quad/']); touch([writePath 'Abs/']);
    quadid = fopen([writePath 'Quad/LeveneVals.tsv'],'w');
    absid = fopen([writePath 'Abs/LeveneVals.tsv'],'w');
    fprintf(quadid,'GroupTested\tSpecimen1\tSpecimen2\tP_value\n'); 
    fprintf(absid,'GroupTested\tSpecimen1\tSpecimen2\tP_value\n');
    %Perform Levene tests if they exist
    for i = 1:length(SpecimenTypes)-1
        for j = (i+1):length(SpecimenTypes)

            quadDict = {}; 
            absDict = {}; 
            quadP = []; absP = [];
            for m = 1:length(totalUnLabels)
                if isempty(dists{i,m}) || isempty(dists{j,m})
                    continue;
                end
                touch([writePath 'Quad/' totalUnLabels{m}]); 
                touch([writePath 'Abs/' totalUnLabels{m}]); 
                if length(dists{i,m}) > length(dists{j,m})
                    dists_i = dists{i,m};
                    dists_j = [dists{j,m};NaN*ones(length(dists{i,m})...
                        -length(dists{j,m}),1)];
                    Specimen1 = SpecimenTypes{i};
                    Specimen2 = SpecimenTypes{j};
                else
                    dists_j = dists{j,m};
                    dists_i = [dists{i,m};NaN*ones(length(dists{j,m})...
                        -length(dists{i,m}),1)];
                    Specimen1 = SpecimenTypes{j};
                    Specimen2 = SpecimenTypes{i};
                end
                P = vartestn([dists_i dists_j],'TestType','LeveneQuadratic',...
                'Display','off');
                save([writePath 'Quad/' totalUnLabels{m} '/' ...
                    SpecimenTypes{i} '_' SpecimenTypes{j} '.mat'],'P');
                quadDict = [quadDict [totalUnLabels{m} '\t' SpecimenTypes{i} '\t'...
                    SpecimenTypes{j} '\t'  num2str(P) '\n']];
                quadP = [quadP P];
                P = vartestn([dists_i dists_j],'TestType','LeveneAbsolute',...
                'Display','off');
                absDict = [absDict [totalUnLabels{m} '\t' SpecimenTypes{i} '\t'...
                    SpecimenTypes{j} '\t'  num2str(P) '\n']];
                absP = [absP P];
                save([writePath 'Abs/' totalUnLabels{m} '/' ...
                    SpecimenTypes{i} '_' SpecimenTypes{j} '.mat'],'P');
            end
            [~,quadInd] = sort(quadP); [~,absInd] = sort(absP);
            for q = quadInd
                fprintf(quadid,quadDict{q});
            end
            for q = absInd
                fprintf(absid,absDict{q});
            end
        end
    end
    fclose(quadid); fclose(absid);
end

%% Intragroup analysis: loop over SpecimenTypes
for i = 1:length(SpecimenTypes)
    load([projectDir SpecimenTypes{i} '/Groups.mat']);
    totalUnLabels = {};
    
    quadDict = {}; absDict = {}; quadP = []; absP = [];
    writePath = [intraStatPath SpecimenTypes{i} '/SampleMean/'];
    %keys{k} '/' curUnLabels{p} '_' curUnLabels{q} ]
    touch(writePath);
    quadPath = [writePath 'Quad/']; touch([writePath 'Quad/']);
    absPath = [writePath 'Abs/']; touch([writePath 'Abs/']);
    quadid = fopen([quadPath 'LeveneVals.tsv'],'w');
    absid = fopen([absPath 'LeveneVals.tsv'],'w');
    fprintf(quadid,'GroupTested\tSubgroup1\tSubgroup2\tP_value\n'); 
    fprintf(absid,'GroupTested\tSubgroup1\tSubgroup2\tP_value\n');
    for k = 1:length(keys)
        touch([quadPath keys{k} '/']); touch([absPath keys{k} '/']);
        currentKey = keys{k};
        totalUnLabels = {};
        curLabels = lower(Groups(currentKey));
        curUnLabels = unique([totalUnLabels;unique(totalLabels{k,i})]);
        dists = cell(1,length(curUnLabels));
        for p = 1:length(curUnLabels)
            inds = find(strcmp(curUnLabels{p},curLabels));
            curMeanVerts = zeros(3,meshCollection{i}{1}.nV);
            for m = inds
                curMesh =meshCollection{i}{m};
                curMeanVerts = curMeanVerts + curMesh.V;
            end
            curMeanVerts = curMeanVerts/length(inds);
            %Project back to Procrustes space
            curMeanVerts = curMeanVerts - repmat(mean(curMeanVerts')',1,size(curMeanVerts,2));
            curMeanVerts = curMeanVerts/norm(curMeanVerts,'fro');
            curDists = zeros(length(inds),1);
            for m = 1:length(inds)
                curDists(m) = norm(meshCollection{i}{inds(m)}.V-curMeanVerts,'fro');
            end
            dists{p} = curDists;
        end
        for p = 1:length(curUnLabels)-1
            for q = p+1:length(curUnLabels)   
                if length(dists{p}) > length(dists{q})
                   dists_p = dists{p};
                   dists_q = [dists{q};NaN*ones(length(dists{p})...
                       -length(dists{q}),1)];
                else
                   dists_q = dists{q};
                   dists_p = [dists{p};NaN*ones(length(dists{q})...
                       -length(dists{p}),1)];
                end
                P = vartestn([dists_p dists_q],'TestType','LeveneQuadratic',...
                'Display','off');
                quadDict = [quadDict [keys{k} '\t' curUnLabels{p} '\t'...
                    curUnLabels{q} '\t'  num2str(P) '\n']];
                quadP = [quadP P];
                save([quadPath 'P.mat'],'P');
                P = vartestn([dists_p dists_q],'TestType','LeveneAbsolute',...
                'Display','off');
                absDict = [absDict [keys{k} '\t' curUnLabels{p} '\t'...
                    curUnLabels{q} '\t'  num2str(P) '\n']];
                absP = [absP P];
                save([absPath 'P.mat'],'P');
            end
        end
        [~,quadInd] = sort(quadP); [~,absInd] = sort(absP);
        for q = quadInd
            fprintf(quadid,quadDict{q});
        end
        for q = absInd
            fprintf(absid,absDict{q});
        end  
    end
    fclose(quadid); fclose(absid);
end
    %Need to extract labels before moving on with analysis so as to
    %make everything consistent
