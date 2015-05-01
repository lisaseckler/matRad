% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad script
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 clear
 close all
 clc

% load patient data, i.e. ct, voi, cst

%load HEAD_AND_NECK
%load TG119.mat
%load PROSTATE.mat
load LIVER.mat
%load BOXPHANTOM.mat

% meta information for treatment plan
pln.SAD             = 1000; %[mm]
pln.isoCenter       = matRad_getIsoCenter(cst,ct,0);
pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [0:72:359]; % [�]
pln.couchAngles     = [0 0 0 0 0]; % [�]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = numel(ct.cube);
pln.voxelDimensions = size(ct.cube);
pln.radiationMode   = 'photons'; % either photons / protons / carbon
pln.bioOptimization = false;   % false indicates physical optimization and true indicates biological optimization
pln.numOfFractions  = 30;

%% initial visualization and change objective function settings if desired
matRadGUI

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);

%% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst,0);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst,0);
end

%% Dose visualization
doseVis = matRad_mxCalcDose(dij,ones(dij.totalNumOfBixels,1),cst);
%matRadGUI

%% inverse planning for imrt
optResult = matRad_fluenceOptimization(dij,cst,pln);
matRadGUI

%% sequencing
if strcmp(pln.radiationMode,'photons')
    %Sequencing = matRad_xiaLeafSequencing(optResult.w,stf,7,1);
    Sequencing = matRad_engelLeafSequencing(optResult.w,stf,7);
    seqResult = matRad_mxCalcDose(dij,Sequencing.w,cst);
    matRad_visCtDose(seqResult,cst,pln,ct);
end

%% dvh and conformity index
matRad_calcDVH(optResult,cst)

%%

% vX = [69 90];
% vY = [69 90];
% vZ = [56 74];
% 
% [i, j, l]=  ind2sub(size(ct.cube),1:1:numel(ct.cube));
%  
% vNew = zeros((vX(2)-vX(1))*(vY(2)-vY(1))*(vZ(2)-vZ(1)),1);
% Counter = 1;
% for cnt = 1:numel(ct.cube)
%     
%  if i(cnt)>=vX(1) && i(cnt)<= vX(2) && ...
%     j(cnt)>=vX(1) && j(cnt)<= vX(2) && ...
%     l(cnt)>=vX(1) && l(cnt)<= vX(2) 
% 
%     vNew(Counter) = cnt;
%     Counter = Counter+1;
% 
%  end
%     
% end
% 
% cst{2,4}=vNew;


% s=whos;
% BytesTot = 0;
% 
% for i = 1:length(s)
%    BytesTot = BytesTot + s(i).bytes;
% end
% load('E:\experimental data\optResultRBExD_liver.mat');
% Slice = optResultRBExD_liver.physicalDose(:,:,120);
% load('E:\experimental data\optResulteffectliver.mat');
% Slice2 = optResulteffectliver.physicalDose(:,:,120);
% 
% diff = abs(Slice-Slice2);
% diff1 = abs(optResultRBExD_liver.RBEWeightedDose(:,:,120)-optResulteffectliver.RBEWeightedDose(:,:,120));
% figure,subplot(221),imshow(optResultRBExD_liver.physicalDose(:,:,120),[]);
%        subplot(222),imshow(diff,[]);
%        subplot(223),imshow(optResultRBExD_liver.RBEWeightedDose(:,:,120),[]);
%        subplot(224),imshow(diff1,[]);



       
       
       
