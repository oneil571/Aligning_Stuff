close all

FPC = {'None','FWER','Bonferroni'};
options.mode = 'native';
radiusType = 'discrete';
%alpha = significance level
%method = type of T-test to use. Values are
%FPC = method of false positive correction
localPath = [statPath 'PairwiseLocal/'];
touch(localPath);
frechMeanDir = 'D://Work/ToothAndClaw/FrechetMean/FrechetMeans/Patch_';
%% Do specified pairwise testing for each specimen type in Groups
% Only makes sense for each indiivudal specimen, cross specimens not
% included
for s = 1:length(SpecimenTypes)
    curPath = [MetaGroupBasePath SpecimenTypes{s} '/'];
    load([curPath 'Groups.mat']);
    load([curPath 'newMeshList.mat']);
    keys = Groups.keys;
    
    %% Loop over the different classification types
    for g = 1:length(Groups.keys)
        curLabels = Groups(keys{g});
        unLabels = unique(curLabels);
        labelGroups = cell(size(unLabels));
        labelGroupInds = cell(size(unLabels));
        %% Separate the meshes into groups
        for i = 1:length(unLabels)
            inds = find(strcmp(curLabels,unLabels{i}));
            labelGroupInds{i} = inds;
            labelGroups{i} = cell(size(inds));
            for j = 1:length(inds)
                labelGroups{i}{j} = newMeshList{inds(j)};
            end
        end
        
        %% Construct mean meshes and average mean mesh
        labelMeans = cell(size(unLabels));
        totalMeanVerts = zeros(size(labelGroups{1}{1}.V));
        for i = 1:length(unLabels)
            len = length(labelGroups{i});
            verts = zeros(size(labelGroups{1}{1}.V));
            for j = 1:length(labelGroups{i})
                verts = verts + labelGroups{i}{j}.V/len;
            end
            labelMeans{i} = Mesh('VF',verts,labelGroups{1}{1}.F);
            labelMeans{i}.V = labelMeans{i}.V-mean(labelMeans{i}.V,2);
            labelMeans{i}.V = labelMeans{i}.V/norm(labelMeans{i}.V,'fro');
            totalMeanVerts = totalMeanVerts + labelMeans{i}.V;
        end
        totalMeanVerts = totalMeanVerts/length(unLabels);
        totalMeanVerts = totalMeanVerts-mean(totalMeanVerts,2);
        totalMeanVerts = totalMeanVerts/norm(totalMeanVerts,'fro');
        totalMean = Mesh('VF',totalMeanVerts,labelGroups{1}{1}.F);
        
        %% Find vertices for each patch based on average mean for classification
        patchInds = cell(totalMean.nV,1);
        disp('Extracting Pathes...')
        switch radiusType
            case 'discrete'
                totalNumWalks = sparse(totalMean.nV,totalMean.nV);
                for k = 1:radius+1
                    totalNumWalks = totalNumWalks + totalMean.A^(k-1);
                end
                for i = 1:totalMean.nV
                    patchInds{i} = find(totalNumWalks(i,:));
                end
            case 'continuous'
                for i = 1:totalMean.nV
                    [D,~,~] = totalMean.PerformFastMarching(i);
                    patchInds{i} = find(D <= raidus);
                end
            otherwise
                error('Invalid radius type input');
        end
        %% Extract Patches based on Vertices and separate by group
        patchList = cell(totalMean.nV,1);
        meanPatches = patchList; euMeans = patchList; adMeans = patchList;
        peroMeans = patchList;
        progressbar
        for i = 1:totalMean.nV
            curPatches = cell(length(newMeshList)+length(unLabels),1);
            triInds = find(sum(ismember(totalMean.F,patchInds{i}))==3);
            [~,faces] = ismember(totalMean.F(:,triInds),patchInds{i});
            for j = 1:length(newMeshList)
                curPatches{j} = Mesh('VF',newMeshList{j}.V(:,patchInds{i}),faces);
                curPatches{j}.V = curPatches{j}.V-mean(curPatches{j}.V,2);
                curPatches{j}.V = curPatches{j}.V/norm(curPatches{j}.V,'fro');
            end
            load([frechMeanDir num2str(i) '/Adapis_Frechet.mat']);
            Frechet_Mean = Frechet_Mean - mean(Frechet_Mean,2);
            Frechet_Mean = Frechet_Mean/norm(Frechet_Mean,'fro');
            adMeans{i} = Mesh('VF',Frechet_Mean,faces);
            curPatches{length(newMeshList)+1} = adMeans{i};
            load([frechMeanDir num2str(i) '/Eulemur_Frechet.mat']);
            Frechet_Mean = Frechet_Mean - mean(Frechet_Mean,2);
            Frechet_Mean = Frechet_Mean/norm(Frechet_Mean,'fro');
            euMeans{i} = Mesh('VF',Frechet_Mean,faces);
            curPatches{length(newMeshList)+2} = euMeans{i};
            load([frechMeanDir num2str(i) '/Perodicticus_Frechet.mat']);
            Frechet_Mean = Frechet_Mean - mean(Frechet_Mean,2);
            Frechet_Mean = Frechet_Mean/norm(Frechet_Mean,'fro');
            peroMeans{i} = Mesh('VF',Frechet_Mean,faces);
            curPatches{length(newMeshList)+3} = peroMeans{i};
            patchList{i} = curPatches;
            meanVerts = (adMeans{i}.V+euMeans{i}.V+peroMeans{i}.V)/3;
            meanPatches{i} = Mesh('VF',meanVerts,faces);
            meanPatches{i}.V = meanPatches{i}.V-mean(meanPatches{i}.V,2);
            meanPatches{i}.V = meanPatches{i}.V/norm(meanPatches{i}.V,'fro');
            progressbar(i/totalMean.nV);
        end
        
        %% Align Patches to Common Template and Extract Distances
        disp('Computing distances');
        dists = cell(totalMean.nV,1);
        switch alignmentMethod
            case 'AverageMean'
                disp('AverageMean template alignment selected, aligning patches');
                progressbar
                for i = 1:totalMean.nV
                    disp(i)
                    curPatches = patchList{i};
                    for j = 1:length(patchList{i})
                        [U,~,V] = svd(curPatches{j}.V*(meanPatches{i}.V'));
                        R = V*(U');
                        newVerts = curPatches{j}.V;
                        for k = 1:size(newVerts,2)
                            newVerts(:,k) = R*newVerts(:,k);
                        end
                        curPatches{j}.V = newVerts;
                    end
                    patchList{i} = curPatches;
                    progressbar(i/totalMean.nV);
                end
                disp('Computing distances')
                progressbar
                for i = 1:totalMean.nV
                    curPatches = patchList{i};
                    curDists = zeros(length(curPatches),length(curPatches));
                    for j = 1:size(curDists,1)
                        for k = 1:size(curDists,2)
                            if j ~=k
                                curDists(j,k) = norm(curPatches{j}.V-curPatches{k}.V,'fro');
                            end
                        end
                    end
                    dists{i} = curDists;
                    progressbar(i/totalMean.nV);
                end
            case 'pairwise'
                disp('Pairwise selected. This method will take a long time');
                disp('Computing distances based on pairwise alignments');
                progressbar
                for i = 1:length(dists)
                    disp(i)
                    curDists = zeros(length(newMeshList),length(newMeshList));
                    for j = 1:length(newMeshList)
                        mesh1 = patchList{i}{j};
                        for k = 1:length(newMeshList)
                            mesh2 = patchList{i}{k};
                            if j ~= k
                                [U,~,V] = svd(mesh1.V*(mesh2.V'));
                                R = V*U';
                                newVerts = mesh1.V;
                                for q = 1:size(newVerts,2)
                                    newVerts(:,q) = R*newVerts(:,q);
                                end
                                mesh1.V = newVerts;
                                curDists(j,k) = norm(mesh1.V-mesh2.V,'fro');
                            end
                        end
                    end
                    dists{i} = curDists;
                    progressbar(i/totalMean.nV);
                end
        end
        
        %% Embed distance functions. Try classical MDS first.
        embedCoords = cell(totalMean.nV,1);
        disp('Computing Embeddings...');
        progressbar
        for i = 1:totalMean.nV
            try
                embedCoords{i} = mdscale(dists{i},embedDim,'Criterion','strain');
            catch
                warning(['Cannot embed patch around vertex ' num2str(i) ' via classical MDS.']);
                try
                    embedCoords{i} = mdscale(dists{i},embedDim);
                catch
                    error('Prescribed distance method not embeddable, please try another metric');
                end
            end
            progressbar(i/totalMean.nV);
        end
        
        %% Setup output directories
        exponents = zeros(size(pValues));
        for p = 1:length(pValues)
            exp = 0;
            while true
                if rem(pValues(p)*10^(exp),1) == 0 
                    break;
                else
                    exp = exp+1;
                end
            end
            exponents(p) = exp;
            %set all folders
            touch([curPath 'HeatMapsMMLS/']);
            touch([curPath 'HeatMapsMMLS/' keys{g} '/']);
            touch([curPath 'HeatMapsMMLS/' keys{g} '/One_Sample/']);
            for c = 1:length(FPC)
                touch([curPath 'HeatMapsMMLS/' keys{g} '/One_Sample/Vertex/' FPC{c} '/']);
                touch([curPath 'HeatMapsMMLS/' keys{g} '/One_Sample/Patch/AverageMean/' alignmentMethod '/' FPC{c} '/']);
                touch([curPath 'HeatMapsMMLS/' keys{g} '/One_Sample/Vertex/' FPC{c} '/'...
                    num2str(pValues(p)*10^(exp)) 'e-' num2str(exp) '/']);
                touch([curPath 'HeatMapsMMLS/' keys{g} '/One_Sample/Patch/AverageMean/' alignmentMethod '/' FPC{c} '/'...
                    num2str(pValues(p)*10^(exp)) 'e-' num2str(exp) '/']);
            end
            for t = 1:length(TTestType)
                touch([curPath 'HeatMapsMMLS/' keys{g} '/' TTestType{t} '/']);
                touch([curPath 'HeatMapsMMLS/' keys{g} '/' TTestType{t} '/Vertex/']);
                touch([curPath 'HeatMapsMMLS/' keys{g} '/' TTestType{t} '/Patch/' alignmentMethod '/']);
                for c = 1:length(FPC)
                    touch([curPath 'HeatMapsMMLS/' keys{g} '/' TTestType{t} '/Vertex/' FPC{c} '/']);
                    touch([curPath 'HeatMapsMMLS/' keys{g} '/' TTestType{t} '/Patch/' alignmentMethod '/Radius_'...
                        num2str(radius) FPC{c} '/']);
                    touch([curPath 'HeatMapsMMLS/' keys{g} '/' TTestType{t} '/Vertex/' FPC{c} '/'...
                        num2str(pValues(p)*10^(exp)) 'e-' num2str(exp) '/']);
                    touch([curPath 'HeatMapsMMLS/' keys{g} '/' TTestType{t} '/Patch/' alignmentMethod '/Radius_' num2str(radius)...
                        '/' FPC{c} '/' num2str(pValues(p)*10^(exp)) 'e-' num2str(exp) '/']);
                    
                end
            end
        end
        
        %% Pairwise Tests
        for i = 1:length(labelGroups)
            for j = (i+1):length(labelGroups)
                %% Set up graphical files
                disp('~~~~~');
                disp([unLabels{i} ' vs ' unLabels{j}]);
                curMeanVerts = .5*(labelMeans{i}.V+labelMeans{j}.V);
                curMeanVerts = curMeanVerts - repmat(mean(curMeanVerts,2),1,labelMeans{i}.nV);
                curMeanVerts = curMeanVerts/norm(curMeanVerts,'fro');
                curMeanMesh = Mesh('VF',curMeanVerts,labelMeans{i}.F);
                groups1 = labelGroups{i}; groups2 = labelGroups{j};
                inds1 = labelGroupInds{i}; inds2 = labelGroupInds{j};
                g1MeanVerts = zeros(size(groups1{1}.V));
                for q = 1: length(groups1)
                    [U,~,V] = svd(groups1{q}.V*(curMeanMesh.V'));
                    R = V*U';
                    newVerts = groups1{q}.V;
                    for s = 1:size(newVerts,2)
                        newVerts(:,s) = R*newVerts(:,s);
                    end
                    groups1{q}.V = newVerts;
                    g1MeanVerts = g1MeanVerts+groups1{q}.V;
                end
                g1MeanVerts = g1MeanVerts/length(groups1);
                g1MeanVerts = g1MeanVerts - repmat(mean(g1MeanVerts,2),1,groups1{1}.nV);
                g1MeanVerts = g1MeanVerts/norm(g1MeanVerts,'fro');
                g1Mean = Mesh('VF',g1MeanVerts,groups1{1}.F);

                g2MeanVerts = zeros(size(groups2{1}.V));
                for q = 1: length(groups2)
                    [U,~,V] = svd(groups2{q}.V*(curMeanMesh.V'));
                    R = V*U';
                    newVerts = groups2{q}.V;
                    for s = 1:size(newVerts,2)
                        newVerts(:,s) = R*newVerts(:,s);
                    end
                    groups2{q}.V = newVerts;
                    g2MeanVerts = g2MeanVerts+groups2{q}.V;
                end
                g2MeanVerts = g2MeanVerts/length(groups2);
                g2MeanVerts = g2MeanVerts - repmat(mean(g2MeanVerts,2),1,groups2{1}.nV);
                g2MeanVerts = g2MeanVerts/norm(g2MeanVerts,'fro');
                g2Mean = Mesh('VF',g2MeanVerts,groups2{1}.F);
                %Perform tests is able
                
                %% Perform 1-sample test and display results
                if length(groups1) == 1 || length(groups2) == 1
                    if max(length(groups1),length(groups2)) <= embedDim+1
                        disp(['Too few samples to do statistical test at desired dimension, skipping']);
                    continue;
                    else
                        disp(['Performing 1-sample T-test'])
                        testPath = [localPath 'OneSample/Patch/' alignmentMethod '/'];
                        touch(testPath);
                        if exist([testPath unLabels{i} '_' unLabels{j} '.mat']) && ~RecomputeLocalStats
                            disp('Raw p-values already computed, reloading...');
                            load([testPath unLabels{i} '_' unLabels{j} '.mat'])
                        else
                            [~,pVal] = TTest_OneSamplePatch(embedCoords,inds1,inds2);
                            save([localPath unLabels{i} '_' unLabels{j} '.mat'],'pVal');
                        end
                        [sortP,idx] = sort(pVal);
                        for p = 1:length(pValues)
                            disp([unLabels{i} ' vs. ' unLabels{j} ' ' TTestType{t} ' ' num2str(pValues(p))]);
                            for c = 1:length(FPC)
                                curFPC = FPC{c};
                                switch curFPC
                                    case 'FWER'
                                        maxP = sortP(max(find(sortP<(pValues(p)*(1:length(sortP))/length(sortP)))));
                                    case 'Bonferroni'
                                        maxP = sortP(max(find(sortP<(pValues(p)/g1Mean.nV))));
                                    case 'None'
                                        maxP = sortP(max(find(sortP<pValues(p))));
                                end
                                if isempty(maxP)
                                    disp(['No significant points detected for ' curFPC ' FPC, skipping']);
                                    continue;
                                end
                                extrPts = double((1-pVal) >= (1-maxP));
                                heatMap = zeros(1,length(extrPts));
                                for k = 1:length(extrPts)
                                    heatMap(idx(k)) = extrPts(idx(k))*(0.2+0.8*(maxP-pVal(idx(k)))/(maxP+0.00001));
                                end
                                clear h
                                clear Link
                                figure
                                h(1) = subplot(1,2,1);
                                g1Mean.ViewFunctionOnMesh(heatMap',options);
                                hold on
                                axis off
                                if viewDisplace
                                    if size(extrPts,1) ~=1
                                        extrPts = extrPts';
                                    end
                                    displace = g2Mean.V - g1Mean.V;
                                    displace(1,:) = displace(1,:).*extrPts;
                                    displace(2,:) = displace(2,:).*extrPts;
                                    displace(3,:) = displace(3,:).*extrPts;
                                    displaceColor = dot(displace,g1Mean.ComputeNormal)/sqrt(sum(displace.^2));
                                    displaceColor = .5*displaceColor+.5;
                                    displaceColorMat = [];
                                    for k = 1:length(displaceColor)
                                        displaceColorMat = [displaceColorMat; [1 1 1]*displaceColor(k)];
                                    end
                                    quiver3(g1Mean.V(1,:),g1Mean.V(2,:),g1Mean.V(3,:),...
                                        displace(1,:),displace(2,:),displace(3,:),'Color','k');
                                end
                                h(2) = subplot(1,2,2);
                                g2Mean.ViewFunctionOnMesh(heatMap',options);
                                hold on
                                axis off
                                if viewDisplace
                                    displace = g1Mean.V - g2Mean.V;
                                    displace(1,:) = displace(1,:).*extrPts;
                                    displace(2,:) = displace(2,:).*extrPts;
                                    displace(3,:) = displace(3,:).*extrPts;
                                    displaceColor = dot(displace,g2Mean.ComputeNormal)/sqrt(sum(displace.^2));
                                    displaceColor = .5*displaceColor+.5;
                                    displaceColorMat = [];
                                    for k = 1:length(displaceColor)
                                        displaceColorMat = [displaceColorMat; [1 1 1]*displaceColor(k)];
                                    end
                                    quiver3(g2Mean.V(1,:),g2Mean.V(2,:),g2Mean.V(3,:),...
                                        displace(1,:),displace(2,:),displace(3,:),'Color','k');
                                end
                                Link = linkprop(h, {'CameraUpVector',...
                                    'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
                                %Quick Conversion To Scientific Notation
                                setappdata(gcf, 'StoreTheLink', Link);
                                figPath = [curPath 'HeatMapsMMLS/' keys{g} '/OneSample/Patch/' alignmentMethod '/' FPC{c} '/' ...
                                    num2str(pValues(p)*10^(exponents(p))) ...
                                    'e-' num2str(exponents(p)) '/'];
                                touch(figPath);
                                saveas(gcf,[figPath unLabels{i} '_' unLabels{j}]);  
                                close all
                            end
                        end
                    end
        %         elseif length(labelGroups{i}) == 1 || length(labelGroups{j}) == o 1
        %             disp(['Too few samples to do family test, must do one sample']);
        %             [~,pVal] = TTest_Standard(labelGroups{i},labelGroups{j});
                %% Skip test if not enough samples
                elseif length(groups1)+length(groups2) <= embedDim+1
                    disp('Not enough samples to do two-sample comparison, skipping');
                
                %% Do pairwise test
                else
                    for t = 1:length(TTestType)
                        disp(TTestType{t});
                        curTest = TTestType{t};
                        testPath = [localPath curTest '/' alignmentMethod '/Radius_' num2str(radius) '/'];
                        touch(testPath)
                        if 0 == 1%exist([testPath unLabels{i} '_' unLabels{j} '.mat']) && ~RecomputeLocalStats
                            disp('Raw p-values already computed, reloading...');
                            load([testPath unLabels{i} '_' unLabels{j} '.mat'])
                        else
                            switch curTest
                                case 'EqualCovariance'
                                    [~,pVal] = TTest_StandardPatchMMLS(embedCoords,inds1,inds2,length(newMeshList)+i,length(newMeshList)+j);
                                case 'UnequalCovariance'
                                    [~,pVal] = TTest_UnequalPatchMMLS(embedCoords,inds1,inds2,length(newMeshList)+i,length(newMeshList)+j);
                                case 'EqualCovariance_Permutation'
                                    [~,pVal] = TTest_Standard_PermutationPatch(embedCoords,inds1,inds2,numPerm);
                                case 'UnequalCovariance_Permutation'
                                    [~,pVal] = TTest_Unequal_PermutationPatch(embedCoords,inds1,inds2,numPerm);
                            end
                            save([testPath unLabels{i} '_' unLabels{j} '.mat'],'pVal');
                        end
                        if size(pVal,1) > 1
                            pVal = pVal';
                        end
                        [sortP,idx] = sort(pVal);
                        for p = 1:length(pValues)
                            disp([unLabels{i} ' vs. ' unLabels{j} ' ' TTestType{t} ' ' num2str(pValues(p))]);
                            for c = 1:length(FPC)
                                curFPC = FPC{c};
                                switch curFPC
                                    case 'FWER'
                                        maxP = sortP(max(find(sortP<(pValues(p)*(1:length(sortP))/length(sortP)))));
                                    case 'Bonferroni'
                                        maxP = sortP(max(find(sortP<(pValues(p)/g1Mean.nV))));
                                    case 'None'
                                        maxP = sortP(max(find(sortP<pValues(p))));
                                end
                                if isempty(maxP)
                                    disp(['No significant points detected for ' curFPC ' FPC, skipping']);
                                    continue;
                                end
                                extrPts = double((1-pVal) >= (1-maxP));
                                heatMap = zeros(1,length(extrPts));
                                for k = 1:length(extrPts)
                                    heatMap(idx(k)) = extrPts(idx(k))*(0.2+0.8*(maxP-pVal(idx(k)))/(maxP+0.00001));
                                end
                                clear h
                                clear Link
                                figure
                                h(1) = subplot(1,2,1);
                                g1Mean.ViewFunctionOnMesh(heatMap',options);
                                hold on
                                axis off
                                if viewDisplace
                                    if size(extrPts,1) ~=1
                                        extrPts = extrPts';
                                    end
                                    displace = g2Mean.V - g1Mean.V;
                                    displace(1,:) = displace(1,:).*extrPts;
                                    displace(2,:) = displace(2,:).*extrPts;
                                    displace(3,:) = displace(3,:).*extrPts;
                                    displaceColor = dot(displace,g1Mean.ComputeNormal)/sqrt(sum(displace.^2));
                                    displaceColor = .5*displaceColor+.5;
                                    displaceColorMat = [];
                                    for k = 1:length(displaceColor)
                                        displaceColorMat = [displaceColorMat; [1 1 1]*displaceColor(k)];
                                    end
                                    quiver3(g1Mean.V(1,:),g1Mean.V(2,:),g1Mean.V(3,:),...
                                        displace(1,:),displace(2,:),displace(3,:),'Color','k');
                                end
                                h(2) = subplot(1,2,2);
                                g2Mean.ViewFunctionOnMesh(heatMap',options);
                                hold on
                                axis off
                                if viewDisplace
                                    displace = g1Mean.V - g2Mean.V;
                                    displace(1,:) = displace(1,:).*extrPts;
                                    displace(2,:) = displace(2,:).*extrPts;
                                    displace(3,:) = displace(3,:).*extrPts;
                                    displaceColor = dot(displace,g2Mean.ComputeNormal)/sqrt(sum(displace.^2));
                                    displaceColor = .5*displaceColor+.5;
                                    displaceColorMat = [];
                                    for k = 1:length(displaceColor)
                                        displaceColorMat = [displaceColorMat; [1 1 1]*displaceColor(k)];
                                    end
                                    quiver3(g2Mean.V(1,:),g2Mean.V(2,:),g2Mean.V(3,:),...
                                        displace(1,:),displace(2,:),displace(3,:),'Color','k');
                                end
                                Link = linkprop(h, {'CameraUpVector',...
                                    'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
                                %Quick Conversion To Scientific Notation
                                setappdata(gcf, 'StoreTheLink', Link);
                                 figPath = [curPath 'HeatMapsMMLS/' keys{g} '/' TTestType{t} '/Patch/' alignmentMethod '/Radius_'...
                                     num2str(radius) '/' FPC{c} '/' num2str(pValues(p)*10^(exponents(p)))...
                                     'e-' num2str(exponents(p)) '/'];
                                touch(figPath);
                                saveas(gcf,[figPath unLabels{i} '_' unLabels{j}]);  
                                close all
                            end
                        end      
                    end
                end
            end
        end
    end
end

