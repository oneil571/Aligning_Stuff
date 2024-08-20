%Makes working directory activated
function touch(workingDir)
    if ~isdir(workingDir)
        mkdir(workingDir);
    end
end