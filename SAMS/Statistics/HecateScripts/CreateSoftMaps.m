%% Load relevant details
load([workingPath 'Names.mat'])
meshList = cell(length(Names),1);
for i = 1:length(Names)
    load([workingPath 'ProcessedMAT/' Names{i} '.mat']);
    meshList{i} = G;
end
touch([workingPath 'SoftMapsMatrix/']);

%% Create full soft maps matrix

progressbar
disp('Creating soft maps matrix');

for i = 1:length(meshList)
    G1 = meshList{i};
    softMapsMatrix = cell(length(meshList),1);
    AugKernel12 = cell(length(meshList),1);
    AugKernel21 = AugKernel12;
    load([workingPath 'TextureCoordsSource/TextureCoordsSource_' num2str(i) '.mat']);
    parfor j = 1:length(meshList)  
        [~,~,AugKernel12{j},~] = MapSoftenKernel(TextureCoordsSource{j}...
            ,meshList{j}.Aux.UniformizationV(1:2,:),...
            meshList{j}.F,G1.V,meshList{j}.V,fiberEps);
        [~,~,AugKernel21{j},~] = MapSoftenKernel(meshList{j}.Aux.UniformizationV(1:2,:)...
            ,TextureCoordsSource{j},G1.F,meshList{j}.V,G1.V,fiberEps);
        softMapsMatrix{j} = max(AugKernel12{j},AugKernel21{j}');
    end
    save([workingPath 'SoftMapsMatrix/SoftMapsMatrix_' num2str(i) '.mat'],'softMapsMatrix');
    clear softMapsMatrix
    clear AugKernel12
    clear AugKernel21
    progressbar(i/length(meshList));
end


