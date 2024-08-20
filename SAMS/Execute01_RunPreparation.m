%CAUTION: The code may still fully run if the meshes in your collection are
%not simply connected. This does not mean that the output is good. No
%guarantees are provided.


%% Initialize and set paths
initialize;
MappingSetup;
disp('Setting Paths');

FORCERECOMPUTATIONOFEVERYTHING = 1;
%{
if FORCERECOMPUTATIONOFEVERYTHING==1 & exists(workingPath)
    rmdir(workingPath) 
end
%}

%% Move data to new directory if not already there.
disp('Initializing Data Collection')
if ~exist([workingPath 'Flags.mat'])
    Flags = containers.Map;
    save([workingPath 'Flags.mat'],'Flags');
else
    load([workingPath 'Flags.mat']);
end
if exist([workingPath 'Names.mat'])
    load ([workingPath 'Names.mat'])
    if length(Names) > 0
        disp('Data already exists in working folder, will only reprocess if forced');
    else
        MoveDataToOutputFolder;
    end
else
    MoveDataToOutputFolder;
end

%% Perform quality check of surfaces
% Code will stop if at least one mesh is not a manifold
% disp('Checking validity of Surfaces');
% if exist([workingPath 'Flags.mat'])
%     load([workingPath 'Flags.mat']);
%     if isKey(Flags,'AreHomeomorphic')
%         if Flags('AreHomeomorphic') == 0
%             rslt = HomeomorphismCheck(workingPath);
%             if rslt.isDisc == -1
%                 Flags('AreHomeomorphic') = 0;
%                 save([workingPath 'Flags.mat'],'Flags');
%                 error('Script stopping, please fix problems with meshes and then try again.');
%             else
%                 Flags('AreHomeomorphic') = 1;
%                 Flags('isDisc') = rslt.isDisc;
%                 save([workingPath 'Flags.mat'],'Flags');
%             end
%         else
%             disp('Already verified homeomorphisms');
%         end
%     else
%         rslt = HomeomorphismCheck(workingPath);
%         if rslt.isDisc == -1
%             Flags('AreHomeomorphic') = 0;
%             save([workingPath 'Flags.mat'],'Flags');
%             error('Script stopping, please fix problems with meshes and then try again.');
%         else
%             Flags('AreHomeomorphic') = 1;
%             Flags('isDisc') = rslt.isDisc;
%             save([workingPath 'Flags.mat'],'Flags');
%         end
%     end
% else
%     rslt = HomeomorphismCheck(workingPath);
%     if rslt.isDisc == -1
%         Flags('AreHomeomorphic') = 0;
%         save([workingPath 'Flags.mat'],'Flags');
%         error('Script stopping, please fix problems with meshes and then try again.');
%     else
%         Flags('AreHomeomorphic') = 1;
%         Flags('isDisc') = rslt.isDisc;
%         save([workingPath 'Flags.mat'],'Flags');
%     end
% end

%% Feature computation for registration step
% Note: errors may come up at this step, unclear why, please report if one
% is found
disp('Computing Necessary Mesh Features');
problemMeshes = zeros(length(Names),1);
if ~isKey(Flags,'FeaturesComputed') || ForceFeatureRecomputation
    if isKey(Flags,'FeaturesComputed')
        disp('Features are already computed, you may safely abort');
    end
    ComputeFeatures;
end
Flags('FeaturesComputed') = 1;
save([workingPath 'Flags.mat'],'Flags');

%% Code successful
disp('Finished preparing meshes, you may begin mapping now');