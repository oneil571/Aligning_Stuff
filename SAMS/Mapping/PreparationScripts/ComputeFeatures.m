%% Set and create paths
touch([workingPath 'ProcessedMAT/']);
output_path = [workingPath 'ProcessedMAT/'];
touch([workingPath 'MappingData/']);
load([workingPath 'RawNames.mat']);
badDataPath = [workingPath '/BadData/'];
badFeatureDir = [badDataPath 'FeatureComputationFailed/'];
touch(badFeatureDir);
load([workingPath 'Flags.mat']);
addpath(genpath('./Mapping/external/AQP_toolbox/'));
%% Get parameters into approrpriate form

%options.isDisc = Flags('isDisc');
options.ConfMaxLocalWidth = ConfMaxLocalWidth;              
options.GaussMaxLocalWidth = GaussMaxLocalWidth;           
options.MeanMinLocalWidth = MeanMinLocalWidth;              
options.DNEMaxLocalWidth = DNEMaxLocalWidth;               
options.SmoothCurvatureFields = SmoothCurvatureFields;                   
options.numGPLmks = numGPLmks;                    
options.pointCloud = 0;

%% Get meshes and do initial, quick feature computations

disp('Loading Meshes...')
%%progressbar
badMeshes = {};
badInds = [];
for i = 1:length(Names)
    if exist([workingPath 'ProcessedMAT/' Names{i} '.mat'],'file') > 0
        continue;
    end
    try
        G = Mesh('off',[workingPath 'RawOFF/' Names{i} '.off']);       %load mesh
        [area,~] = G.ComputeSurfaceArea; G.V = G.V/sqrt(real(area));
        [~,TriAreas] = G.ComputeSurfaceArea;
        while true
            delInds = [];
            if max(abs(imag(TriAreas))) > 0
                delInds = [delInds;find(abs(imag(TriAreas))>0)];
            end
            if min(TriAreas)<nullFaceBound
                delInds = [delInds;find(TriAreas < nullFaceBound)];
            end
            if isempty(delInds)
                break;
            else
                
                G.F(:,delInds) = [];
                G = Mesh('VF',G.V,G.F);
                delInds = find(~ismember(1:G.nV,reshape(G.F,1,3*G.nF)));
                G.DeleteVertex(delInds);
                [bd,~] = G.FindOrientedBoundaries;
                if min(size(bd)>1)
                    badMeshes = [badMeshes Names{i}];
                    badInds = [badInds i];
                    break;
                end
                [~,TriAreas] = G.ComputeSurfaceArea;
            end
        end

        G.Write([workingPath 'RawOFF/' Names{i} '.off'],'off',options);
        save([workingPath 'RawMAT/' Names{i} '.mat'],'G');
        G.Nf = G.ComputeFaceNormals; %Normal at a face
        G.Nv = G.F2V'*G.Nf';         %Defining normal vector at a vertex
        G.Nv = G.Nv'*diag(1./sqrt(sum((G.Nv').^2,1)));
        G.nF = size(G.F,2);
        G.nV = size(G.V,2);     %needed in case deletion after check
        %%progressbar(i/length(Names));
    catch
        badInds = [badInds i];
        %progressbar(i/length(Names));
    end
end


%% Do complex feature computations and save
disp('Computing surface features');
%%progressbar
problemMeshes = zeros(length(Names),1);
badMeshList = {};
badMeshInds = [];

%% Loop over meshes to extract features
numDiscs = 0;
options.isDisc = 0;
for i = 1:length(Names)%1:length(Names)
    if ismember(i,badInds)
        disp(['Error in mesh ' num2str(i)]);
        badMeshInds = [badMeshInds i];
        badMeshList = [badMeshList Names{i}];
        continue
    end
    if exist([workingPath 'ProcessedMAT/' Names{i} '.mat'],'file') > 0
        disp('Mesh already processed, continuing...')
        continue;
    end
    G = Mesh('off',[workingPath 'RawOFF/' Names{i} '.off']);
    %%progressbar(i/length(Names));
    % this is dumb but will work for now
     

    try
        curBoundary = G.FindOrientedBoundaries();
        if min(size(curBoundary)) > 1
            disp('Current mesh is not simply connected, skipping');
            error();
        elseif min(size(curBoundary)) == 1
            numDiscs = numDiscs+1;
            options.isDisc = 1;
        else
            options.isDisc = 0;
        end
        [G,Aux] =G.ComputeAuxiliaryInformation(options);
        G.Aux = Aux;
        G.Aux.Name = Names{i};
        save([output_path Names{i} '.mat'],'G');
    catch
        disp(['Error in mesh ' num2str(i)]);
        badMeshInds = [badMeshInds i];
        badMeshList = [badMeshList Names{i}];
    end
end

%% Repeat computations

if ~isempty(badMeshList)
    disp('ALERT: The following meshes require manual cleaning to continue')
    for i = 1:length(badMeshList)
        disp(badMeshList{i})
        copyfile([workingPath 'RawOFF/' badMeshList{i} '.off'],[badFeatureDir badMeshList{i} '.off']);
    end
    Names(badMeshInds) = [];
    disp('These meshes will not be used');
    disp('Please edit the files in the RawOFF output directory to use them');
end
if numDiscs == 0  && isempty(badMeshInds)
    Flags('isDisc') = 0;
elseif numDiscs == length(Names) && isempty(badMeshInds)
    Flags('isDisc') = 1;
else
    disp('Not all meshes have the same topology. Please check your data to verify whether you want disk or sphere topology.');
end
save([workingPath 'Names.mat'],'Names');
save([workingPath 'Flags.mat'],'Flags');

%% Saving

