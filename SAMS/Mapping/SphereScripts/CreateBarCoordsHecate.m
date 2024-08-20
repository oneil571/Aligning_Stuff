%% Setup and loading
rmpath(genpath([SAMSPath '/utils/']));
addpath(genpath([SAMSPath 'Mapping/external/']));
DataDir = [workingPath 'OrbifoldDataHecate/'];
NamesPath = [workingPath 'Names.mat'];
load(NamesPath);

%% Flatten meshes
flatteners = {};
for i = 1:length(Names)
    [V,T]=read_off([DataDir Names{i} '.off']);
    inds=load([DataDir Names{i} '.txt']);
    flattener1=Flattener(V,T,inds);
    if i == 1
        flattener1.orderTS();
        flattener1.flatten_orbifold();
        %fix numerical errors if exist
        flattener1.fixFlipsNew();
        %add to the cell array of flatteners
        leadFlattener=flattener1;
    else
        flattener1.uncut_cone_inds=flattener1.uncut_cone_inds(leadFlattener.reorder_cones);
        flattener1.flatten_orbifold();
        %fix numerical errors if exist
        flattener1.fixFlipsNew();
            %add to the cell array of flatteners
    end
    flatteners{end+1}=flattener1;
end

%% Compute maps
map=UncutSurfMap(flatteners);
disp('Computing Maps')
progressbar
for i = 1:length(Names)
    for j = 1:length(Names)
        if i ~= j
            map.compute(i,j)
        end
        progressbar(((i-1)*length(Names)+j)/(length(Names)^2));
    end
end
map = map.barCoords;
for i = 1:length(Names)
    map{i,i} = sparse(eye(size(map{1,2},1)));
end
save([workingPath 'HecateMaps.mat'],'map');

