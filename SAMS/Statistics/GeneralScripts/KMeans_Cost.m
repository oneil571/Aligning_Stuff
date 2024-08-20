function val = KMeans_Cost(data, classA)
% data - column vecs
% ClassA - indices for class A in TF values
avgA = mean(data(:, classA), 2);
avgB = mean(data(:, setdiff(1:size(data,2),classA)), 2);
sumA = sum(vecnorm(data(:,classA) - repmat(avgA, 1, length(classA)), 2).^2);
sumB = sum(vecnorm(data(:,setdiff(1:size(data,2),classA))...
    - repmat(avgB, 1, size(data,2)-length(classA)), 2).^2);
val = sumA+sumB;
end