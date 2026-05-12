%% Simplifying HIT base data
% only taking the third focus Index to make sure the TopasMCEngine has no
% problem in line 1692

load("carbon_HITfluenceSpectra_updatedBig.mat")

for i = 1:255
    machine.data(i).initFocus.dist = machine.data(i).initFocus.dist(3,:);
    machine.data(i).initFocus.sigma = machine.data(i).initFocus.sigma(3,:);
    machine.data(i).initFocus.emittance = machine.data(i).initFocus.emittance(3);
end

machine.meta.machine = "carbon_HITfluenceSpectra_only1FocusIx";

%% saving
save("C:\Users\l813r\Documents\GitHub\LisaSeckler\matRad\matRad\basedata\carbon_HITfluenceSpectra_only1FocusIx.mat","machine")