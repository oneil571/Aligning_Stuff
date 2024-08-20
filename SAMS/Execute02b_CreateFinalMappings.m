disp('Interpolating sparse correspondences...')
load([workingPath 'Flags.mat']);
if Flags('isDisc') == 0
    SetupHypOrb;
    CreateFinalMappingsSphere;
else
    CreateFinalMappingsDisc2;
end

disp('Mappings computed. Please visualize with plotColorMap before continuing');