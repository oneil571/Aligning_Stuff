%% Load data and extract set number of GP Landmarks for matching
namesPath = [workingPath 'Names.mat'];
MATPath = [workingPath 'ProcessedMAT/'];
mappingPath = [workingPath 'MappingData/'];
load(namesPath);


meshList = cell(length(Names),1);
GPLmkList = cell(length(Names),1);
GPPtClouds = cell(length(Names),1);
for i = 1:length(Names)
    load([MATPath Names{i} '.mat']);
    meshList{i} = G;
    GPLmkList{i} = G.Aux.GPLmkInds;
    GPPtClouds{i} = G.V;
    %centralize point clouds
    GPPtClouds{i} = GPPtClouds{i} - repmat(mean(GPPtClouds{i}'),size(GPPtClouds{i},2),1)';
    GPPtClouds{i} = GPPtClouds{i}/norm(GPPtClouds{i});
end

%% First step: load distances if they exist or compute if needed


if Flags('hasDists')
    disp('Distances already found, loading...');
    load([workingPath 'GPDists.mat']);
    disp('Distances loaded');
else
    GPDists = zeros(length(Names),length(Names));
    disp('Distances not found, computing...')
    progressbar
    for i = 1:length(Names)
        procMaps = cell(length(Names),length(Names));
        dummy = cell(length(Names),1);
        curCloud = GPPtClouds{i};
        parfor j = 1:length(Names)
            if i ~=j
                %[P,procDists(i,j),~] = linassign(ones(size(GPPtClouds{i},2),size(GPPtClouds{i},2)),D);
                %Get map from permutation
                [~,dummy{j}] = knnsearch(GPPtClouds{j}',...
                    curCloud');
                GPDists(i,j) = mean(dummy{j});
            end
        end
        progressbar(i/length(Names));
    end
    GPDists = (GPDists+GPDists')/2;
    save([workingPath 'GPDists.mat'],'GPDists');
    Flags('hasDists') = 1;
    save([workingPath 'Flags.mat','Flags']);
end

if ~isKey(Flags,'hasFlows')
    frInd = find(sum(GPDists.^2)==min(sum(GPDists.^2)));
    frInd = frInd(1);
    Flows = ComputeDirectedFlows(GPDists,1:size(GPDists,1),frInd);
    Flags('hasFlows') = 1;
    save([workingPath 'Flags.mat'],'Flags');
    save([workingPath 'Flows.mat'],'Flows');
else
    frInd = find(sum(GPDists.^2)==min(sum(GPDists.^2)));        %This is needed later
    disp('Flows to Frechet Mean already computed, loading...');
    load([workingPath 'Flows.mat']);
end

            