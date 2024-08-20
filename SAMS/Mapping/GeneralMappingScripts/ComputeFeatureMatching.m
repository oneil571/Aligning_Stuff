meshList = cell(length(Names),1);
for i = 1:length(Names)
    load([workingPath 'ProcessedMAT/' Names{i} '.mat']);
    meshList{i} = G;
end

load([workingPath 'GPDists.mat']);
frechMean = find(min(sum(GPDists.^2))==sum(GPDists.^2));

featureList = cell(length(Names),1);
switch featureMap
    case 'Conf'
        for i = 1:length(Names)
            featureList{i} = meshList{i}.Aux.ConfMaxInds;
        end
    case 'Gauss'
        for i = 1:length(Names)
            featureList{i} = meshList{i}.Aux.GaussMaxInds;
        end
    case 'Mean'
        for i = 1:length(Names)
            featureList{i} = meshList{i}.Aux.MeanMinInds;
        end
    case 'DNE'
        for i = 1:length(Names)
            featureList{i} = meshList{i}.Aux.DNEMaxInds;
        end
end
featureMatchesPairs = cell(length(Names),1);
for i = 1:length(Names)
    disp(['Computing features for ' Names{i}]);
    if i ~= frechMean
        pairMatches = [];
        map_12 = knnsearch(meshList{frechMean}.V(:,featureList{frechMean})',...
            meshList{i}.V(:,featureList{i})');
        map_21 = knnsearch(meshList{i}.V(:,featureList{i})',...
            meshList{frechMean}.V(:,featureList{frechMean})');
        for j = 1:length(map_12)
            
            %You did this because you wanted to restrict analysis to GP
            %landmarks, so you projected onto GP landmarks because you were
            %concerned that there would be overlap
            if map_21(map_12(j)) == j
                lmk12 = knnsearch(meshList{frechMean}.V',...
            meshList{i}.V(:,featureList{i}(j))');
                lmk21 = knnsearch(meshList{i}.V',...
            meshList{frechMean}.V(:,featureList{frechMean}(map_12(j)))');
                weightAdj1 = pdist2(meshList{i}.V',meshList{i}.V').*meshList{i}.A;
                weightAdj2 = pdist2(meshList{frechMean}.V',meshList{frechMean}.V')...
                    .*meshList{frechMean}.A;
                [dist12,~,~] = graphshortestpath(weightAdj2,lmk12,...
                    featureList{frechMean}(map_12(j)));
                [dist21,~,~] = graphshortestpath(weightAdj1,lmk21,...
                   featureList{i}(j));
                
                if dist12 < maxDistTol && dist21 < maxDistTol
                    featureMatchesPairs{i} = [featureMatchesPairs{i};...
                        featureList{i}(j) featureList{frechMean}(map_12(j))];
                  
                end
            end
        end
    end
end
Flags('featureMappings') = 1;
save([workingPath 'Flags.mat'],'Flags');
save([workingPath 'MappingData/FeatureMatches.mat'],'featureMatchesPairs');