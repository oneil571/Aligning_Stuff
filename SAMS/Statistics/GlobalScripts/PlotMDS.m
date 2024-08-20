function PlotMDS(workingPath,plotEllipse,SkipLabels)
close all
if nargin > 2
    testLabels = 1;
else
    testLabels = 0;
    if nargin < 2
        plotEllipse = 0;
    end
end
load([workingPath 'Groups.mat']);
load([workingPath 'MDSEmbedding.mat']);
%Y = mdscale(distMatrix([1:20 22:49],[1:20 22:49]),2);
keys = Groups.keys;
for g = 1:length(keys)
    close all
    labels = Groups(keys{g});
    unLabels = unique(labels);
    colors = distinguishable_colors(length(unLabels));

    figure; hold on;
    for i = 1:length(unLabels)
        inds = find(strcmp(labels,unLabels{i}));
        scatter3(Y(inds,1),Y(inds,2),Y(inds,3),100,colors(i,:),'filled');

    end
    if plotEllipse
        for i = 1:length(unLabels)
            if testLabels
                skipFlag = 0;
                for ll = 1:length(SkipLabels)
                    if strcmp(unLabels{i},SkipLabels{ll})
                        skipFlag = 1;
                        break;
                    end
                end
                if skipFlag
                    continue;
                end
            end
            inds = find(strcmp(labels,unLabels{i}));
            if(length(inds)>2)
                [coeffs,~,latent] = pca(Y(inds,:));
                if(length(latent) <3)
                    coeffs = [coeffs, cross(coeffs(:,1),coeffs(:,2))];
                    [x,y,z] = ellipsoid(0,0,0,sqrt(latent(1)),sqrt(latent(2)),eps);
                else
                    [x,y,z] = ellipsoid(0,0,0,sqrt(latent(1)),sqrt(latent(2)),sqrt(latent(3)));
                end
                newX = reshape(x,1,21^2);
                newY = reshape(y,1,21^2);
                newZ = reshape(z,1,21^2);
                coords = [newX;newY;newZ];
                for j = 1:size(coords,2)
                    coords(:,j) = coeffs*coords(:,j);
                end
                center = mean(Y(inds,:));
                coords = coords+repmat(center',1,21^2);
                x = reshape(coords(1,:),21,21);
                y = reshape(coords(2,:),21,21);
                z = reshape(coords(3,:),21,21);
                surf(x,y,z,'EdgeColor','none','FaceColor',colors(i,:),'FaceAlpha',.5);
            else
                plot3(Y(inds,1),Y(inds,2),Y(inds,3),'Color',colors(i,:));
            end
        end
    end
legend(unLabels);
touch([workingPath '/MDS/']);
figName = ['MDS_Group_' keys{g} '_'];
if plotEllipse
    figName = [figName 'WithEllipse_'];
else
    figName = [figName 'NoEllipse_'];
end

if testLabels
    SkipStr = 'Skip_';
    for ll=1:length(SkipLabels)
        SkipStr = [SkipStr SkipLabels{ll}];
        if ll < length(SkipLabels)
            SkipStr = [SkipStr '_'];
        end
    end
    figName = [figName SkipStr '/'];
else
    figName = [figName 'NoSkip/'];
end
touch([workingPath '/MDS/' figName]);
figDir = [workingPath '/MDS/' figName];
savefig([figDir 'MDS_Normal.fig']);
saveas(gcf,[figDir 'MDS_Normal.png']);
end
