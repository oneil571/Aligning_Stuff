function result = HomeomorphismCheck(workingPath)

%This function serves two purposes. It will first check each surface in a
%particular folder to make sure that it is a manifold. It will
%simultaneously check to see if the surfaces are or are not discs.
%Nonmanifold surfaces and topology discrepancies will be reported.

%Function will return 1 if all manifolds are discs, 0 if all manifolds are
%not discs, and -1 if discrepancies are reported.
badDataPath = [workingPath '/BadData/'];
nonManifoldMeshes = [badDataPath 'NonManifoldMeshes/'];
touch(nonManifoldMeshes);
load([workingPath 'RawNames.mat']);
problemMeshes = {};
badBoundaryMeshes = {};
discBoundaryMeshes = {};
nonDiscBoundaryMeshes = {};
numDiscs = 0; numNonDiscs = 0;
isMan = 2;
progressbar
for i = 1:length(Names)
    disp(Names{i})
    G = Mesh('off',[workingPath 'RawOFF/' Names{i} '.off']);
    delInds = [];
    for q = 1:G.nV
        if length(find(sum(G.F==q))) == 0
            delInds = [delInds q];
        end
    end
    G.DeleteVertex(delInds);
    save([workingPath 'RawMAT/' Names{i} '.mat'],'G');
    
    %Verify that there are no degenerate faces
    badFaceList = [];
    for q = 1:G.nF
        if length(G.F(:,q)) ~= length(unique(G.F(:,q)))
            badFaceList = [badFaceList q];
        end
    end
    G.F(:,badFaceList) = [];
    G = Mesh('VF',G.V,G.F);
    save([workingPath 'RawMAT/' Names{i} '.mat'],'G');

    try
        isManifoldResult = isManifold(G);
        G.FindOrientedBoundaries;
    catch
        isManifoldResult.manifold = 0;
    end
    %Check if manifold. If not, add to list. If yes, check boundary
    if ~(isManifoldResult.manifold == 1)
        problemMeshes = [problemMeshes Names{i}];
    else
        [boundary,~] = G.FindOrientedBoundaries;
        %isManifoldResult.manifold = 0;
        if min(size(boundary)) == 0
            numNonDiscs = numNonDiscs+1;
            nonDiscBoundaryMeshes = [nonDiscBoundaryMeshes Names{i}];
        elseif min(size(boundary)) == 1
            numDiscs = numDiscs+1;
            discBoundaryMeshes = [discBoundaryMeshes Names{i}];
        else
            problemMeshes = [problemMeshes Names{i}];
            badBoundaryMeshes = [badBoundaryMeshes Names{i}];
        end
    end
    progressbar(i/(length(Names)))
end

result.problemMeshes = problemMeshes;
result.badBoundaryMeshes = badBoundaryMeshes;
result.discTopologyMeshes = discBoundaryMeshes;
result.nonDiscTopologyMeshes = nonDiscBoundaryMeshes;
result.numNonDiscs = numNonDiscs;
result.numDiscs = numDiscs;

%In future this should be replaced with code that checks homology
if ~isempty(badBoundaryMeshes)
    disp('ALERT: The manifold meshes do not have consistent topology. Please resolve before proceeding.');
    disp(['There are ' num2str(numDiscs) ' of those meshes with disc topology and ' ...
        num2str(numNonDiscs) ' with non-disc topology.']);
    disp(['There are ' num2str(length(badBoundaryMeshes)) ...
        ' meshes with boundary problems. Please check appropriately.']);
    result.isDisc = -1;
    isMan = -1;
end

if ~isempty(problemMeshes)
    disp('ALERT: The following meshes must be cleaned before proceeding:')
    for i = 1:length(result.problemMeshes)
        disp(result.problemMeshes{i});
    end
    result.isDisc = -1;
    isMan = -1;
end

%% Print results to file
fid = fopen([workingPath 'problemMeshes.txt'],'w');
for i = 1:length(problemMeshes)
    fprintf(fid,[result.problemMeshes{i} '\n']);
end
fclose(fid);

fid = fopen([workingPath 'discTopologyMeshes.txt'],'w');
for i = 1:length(result.discTopologyMeshes)
    fprintf(fid,[result.discTopologyMeshes{i} '\n']);
end
fclose(fid);

fid = fopen([workingPath 'nonDiscTopologyMeshes.txt'],'w');
for i = 1:length(result.nonDiscTopologyMeshes)
    fprintf(fid,[result.nonDiscTopologyMeshes{i} '\n']);
end
fclose(fid);

%% Remaining output
if ~(isMan == -1)
    if numNonDiscs > 0
        result.isDisc = 0;
    else
        result.isDisc = 1;
    end
end

end