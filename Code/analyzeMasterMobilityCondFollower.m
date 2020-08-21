function [fracTracksWithF,fracTracksWithFperMode,fracTracksDiffModeCombo,...
    fracTracksDiffModeNoF,meanDiffCoefBeforeModeCombo,meanDiffCoefDuringModeCombo,...
    ratioDiffCoefDurBefModeCombo,meanDiffCoefModeNoF] = ...
    analyzeMasterMobilityCondFollower(tracksMF,tracksF,diffModeDividerStruct)
%ANALYZEMASTERMOBILITYCONDFOLLOWER analyzes mobility of master channel molecules conditional on follower channel properties
%
%SYNPOSIS [fracTracksWithF,fracTracksWithFperMode,fracTracksDiffModeCombo,...
%    fracTracksDiffModeNoF,meanDiffCoefBeforeModeCombo,meanDiffCoefDuringModeCombo,...
%    ratioDiffCoefDurBefModeCombo,meanDiffCoefModeNoF] = ...
%    analyzeMasterMobilityCondFollower(tracksMF,tracksF,diffModeDividerStruct)
%
%INPUT  tracksMF, tracksF    : Output of getFollowerTracksFromMaster.
%                              MUST BE IN ALTERNATIVE FORMAT (see
%                              getFollowerTracksFromMaster for format details).
%       diffModeDividerStruct: Diffusion mode divider structure, as input
%                              to trackDiffModeAnalysis.
%
%OUTPUT fracTracksWithF             : Fraction of tracks with significant
%                                     follower presence among all tracks.
%       fracTracksWithFperMode      : Fraction of tracks with significant
%                                     follower presence per diffusion mode
%                                     (among all tracks with that diffusion
%                                     mode). 
%                                     Modes are listed from fastest to
%                                     slowest, here and for all variables
%                                     involving information per mode.
%       fracTracksDiffModeCombo     : (Number of modes)-by-(number of
%                                     modes) matrix indicating fraction of
%                                     tracks with various mode
%                                     combinations for before (row) and
%                                     during (column) follower presence,
%                                     among tracks with follower.
%       fracTracksDiffModeNoF       : Fraction of tracks per mode among
%                                     tracks without a follower.
%       meanDiffCoefBeforeModeCombo : Matrix storing mean diffusion
%                                     coefficient before follower
%                                     presence.
%       meanDiffCoefDuringModeCombo : Matrix storing mean diffusion
%                                     coefficient during follower
%                                     presence.
%       ratioDiffCoefDurBefModeCombo: Matrix storing ratio of
%                                     diffusion coefficient during to
%                                     before follower presence.
%       meanDiffCoefModeNoF         : Vector of mean diffusion coefficient
%                                     per mode for tracks without follower
%                                     presence.
%
%Khuloud Jaqaman, July 2019
%Zachariah Malik, May 2020, account for case of no followers
%
% Copyright (C) 2020, Jaqaman Lab - UT Southwestern 
%
% This file is part of timeCourseAnalysisPlus.
% 
% timeCourseAnalysisPlus is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% timeCourseAnalysisPlus is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with timeCourseAnalysisPlus.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

%% Track separation based on follower presence

%convert to matrix format for ease of handling
tracksMat = convStruct2MatIgnoreMS(tracksMF);
[~,followerIndxMat] = convStruct2MatIgnoreMS(tracksF);
followerIndxMat = (followerIndxMat ~= 0); %logical indicating whether follower exists (1) or not (0)

%separate master tracks based on follower properties
[tracksMatNoF,tracksMatIntermF,tracksMatWithFBef,~,tracksMatWithFDuring] = ...
    separateMasterTracksBasedOnFollower(tracksMat,followerIndxMat);

%fraction of tracks that have significant follower presence (among all
%tracks)
% if fracTracksWithF = 0, then no tracks currently have a follower
fracTracksWithF = size(tracksMatWithFDuring,1)/(size(tracksMat,1));

%%%ZMALIK 20200528
% Check if tracks have a follower
if isempty(tracksMatIntermF(~isnan(tracksMatIntermF)))
  anyTrackIntermWithF = 0;
else
  anyTrackIntermWithF = 1;
end
if isempty(tracksMatWithFBef)
  anyTrackWithFBef = 0;
else
  anyTrackWithFBef = 1;
end
%%%
%% Diffusion mode analysis

%%%ZMALIK 20200528 Update following variables depending on whether or not
%%%current set has followers.
% First consider special case of tracks not having a follower
if anyTrackIntermWithF == 0
%    tracksIntermF = struct('tracksFeatIndxCG',0,'tracksCoordAmpCG',NaN,'seqOfEvents',NaN);
    diffModeResIntermF = repmat(struct('diffMode',NaN,'diffCoef',NaN,'msdF2F',NaN,...
       'meanPosStd',NaN,'diffRadius',NaN,'lifetime',NaN),size(tracksMat,1),1);
    diffModeIntermF = vertcat(diffModeResIntermF.diffMode);
else
    tracksIntermF = convertMat2Struct2(tracksMatIntermF);
    diffModeResIntermF = trackDiffModeAnalysis(tracksIntermF,diffModeDividerStruct);
    diffModeIntermF = vertcat(diffModeResIntermF.diffMode);
end

if anyTrackWithFBef == 0
%    tracksWithFBeg = struct('tracksFeatIndxCG',0,'tracksCoordAmpCG',NaN,'seqOfEvents',NaN);
    diffModeResWithFBef = repmat(struct('diffMode',NaN,'diffCoef',NaN,'msdF2F',NaN,...
       'meanPosStd',NaN,'diffRadius',NaN,'lifetime',NaN),size(tracksMat,1),1);
    diffModeWithFBef = vertcat(diffModeResWithFBef.diffMode);
    diffCoefWithFBef = vertcat(diffModeResWithFBef.diffCoef);
else
    tracksWithFBef = convertMat2Struct2(tracksMatWithFBef);
    diffModeResWithFBef = trackDiffModeAnalysis(tracksWithFBef,diffModeDividerStruct);
    diffModeWithFBef = vertcat(diffModeResWithFBef.diffMode);
    diffCoefWithFBef = vertcat(diffModeResWithFBef.diffCoef);
end

if fracTracksWithF == 0
%    tracksWithFDuring = struct('tracksFeatIndxCG',0,'tracksCoordAmpCG',NaN,'seqOfEvents',NaN);
    diffModeResWithFDuring = repmat(struct('diffMode',NaN,'diffCoef',NaN,'msdF2F',NaN,...
       'meanPosStd',NaN,'diffRadius',NaN,'lifetime',NaN),size(tracksMat,1),1);
    diffModeWithFDuring = vertcat(diffModeResWithFDuring.diffMode);
    diffCoefWithFDuring = vertcat(diffModeResWithFDuring.diffCoef);
else
    tracksWithFDuring = convertMat2Struct2(tracksMatWithFDuring);
    diffModeResWithFDuring = trackDiffModeAnalysis(tracksWithFDuring,diffModeDividerStruct);
    diffModeWithFDuring = vertcat(diffModeResWithFDuring.diffMode);
    diffCoefWithFDuring = vertcat(diffModeResWithFDuring.diffCoef);
end

% Now deal with tracks with no follower presence
tracksNoF = convertMat2Struct2(tracksMatNoF);
diffModeResNoF = trackDiffModeAnalysis(tracksNoF,diffModeDividerStruct);
diffModeNoF = vertcat(diffModeResNoF.diffMode);
diffCoefNoF = vertcat(diffModeResNoF.diffCoef);

%%%

% %convert from matrices to structures to call diffusion mode analysis
% tracksNoF = convertMat2Struct2(tracksMatNoF);
% tracksIntermF = convertMat2Struct2(tracksMatIntermF);
% tracksWithFBef = convertMat2Struct2(tracksMatWithFBef);
% % tracksWithFAft = convertMat2Struct2(tracksMatWithFAft);
% tracksWithFDuring = convertMat2Struct2(tracksMatWithFDuring);
% % tracksWithFBetween = convertMat2Struct2(tracksMatWithFBetween);
% 
% %call diffusion mode analysis
% diffModeResNoF = trackDiffModeAnalysis(tracksNoF,diffModeDividerStruct);
% diffModeResIntermF = trackDiffModeAnalysis(tracksIntermF,diffModeDividerStruct);
% diffModeResWithFBef = trackDiffModeAnalysis(tracksWithFBef,diffModeDividerStruct);
% % diffModeResWithFAft = trackDiffModeAnalysis(tracksWithFAft,diffModeDividerStruct);
% diffModeResWithFDuring = trackDiffModeAnalysis(tracksWithFDuring,diffModeDividerStruct);
% % diffModeResWithFBetween = trackDiffModeAnalysis(tracksWithFBetween,diffModeDividerStruct);
% 
% %assemble diffusion mode information per track
% diffModeNoF = vertcat(diffModeResNoF.diffMode);
% diffModeIntermF = vertcat(diffModeResIntermF.diffMode);
% diffModeWithFBef = vertcat(diffModeResWithFBef.diffMode);
% % diffModeWithFAft = vertcat(diffModeResWithFAft.diffMode);
% diffModeWithFDuring = vertcat(diffModeResWithFDuring.diffMode);
% % diffModeWithFBetween = vertcat(diffModeResWithFBetween.diffMode);
% 
% %assemble diffusion coefficient information per track
% diffCoefNoF = vertcat(diffModeResNoF.diffCoef);
% % diffCoefIntermF = vertcat(diffModeResIntermF.diffCoef);
% diffCoefWithFBef = vertcat(diffModeResWithFBef.diffCoef);
% % diffCoefWithFAft = vertcat(diffModeResWithFAft.diffCoef);
% diffCoefWithFDuring = vertcat(diffModeResWithFDuring.diffCoef);
% % diffCoefWithFBetween = vertcat(diffModeResWithFBetween.diffCoef);

%calculate fraction of tracks with follower presence based on motion mode
%before follower (among all tracks)
%this expands on the variable fracTracksWithF calculated above
numMode = size(diffModeDividerStruct.divider,3) + 1;
fracTracksWithFperMode = NaN(numMode,1);
diffModeNoIntermF = [diffModeNoF; diffModeIntermF];
for iMode = 1 : numMode
    indxWithF = find( diffModeWithFBef==(numMode-iMode+1) );
    indxNoIntermF = find( diffModeNoIntermF==(numMode-iMode+1) );
    fracTracksWithFperMode(iMode) = length(indxWithF) / (length(indxNoIntermF)+length(indxWithF));
end

%identify tracks with various mode combinations of before and during
%follower presence
%store number of tracks with these combinations
indxModeCombo = cell(numMode);
numDiffModeCombo = zeros(numMode);
for iMode = 1 : numMode
    for jMode = 1 : numMode
        indxTmp = find( diffModeWithFBef==(numMode-iMode+1) & diffModeWithFDuring==(numMode-jMode+1) );
        indxModeCombo{iMode,jMode} = indxTmp;
        numDiffModeCombo(iMode,jMode) = length(indxTmp);
    end
end

%get fraction of tracks with these various mode combinations (among tracks
%with follower)
fracTracksDiffModeCombo = numDiffModeCombo / sum(numDiffModeCombo(:));

%get diffusion coefficient during and diffusion coefficient before,
%for the different mode combinations
%also get the during-to-before ratio
[meanDiffCoefBeforeModeCombo,meanDiffCoefDuringModeCombo,ratioDiffCoefDurBefModeCombo] = deal(zeros(numMode));
for iMode = 1 : numMode
    for jMode = 1 : numMode
        diffCoefBef = diffCoefWithFBef(indxModeCombo{iMode,jMode});
        diffCoefDur = diffCoefWithFDuring(indxModeCombo{iMode,jMode});
        meanDiffCoefBeforeModeCombo(iMode,jMode) = mean(diffCoefBef);
        meanDiffCoefDuringModeCombo(iMode,jMode) = mean(diffCoefDur);
        ratioDiffCoefDurBefModeCombo(iMode,jMode) = mean(diffCoefDur./diffCoefBef);
    end
end

%get fraction of tracks per mode among tracks without follower
%also get mean diffusion coefficient per mode for tracks without follower
[numDiffModeNoF,meanDiffCoefModeNoF] = deal(zeros(numMode,1));
for iMode = 1 : numMode
    indxModeNoF = find( diffModeNoF==(numMode-iMode+1) );
    numDiffModeNoF(iMode) = length(indxModeNoF);
    meanDiffCoefModeNoF(iMode) = mean(diffCoefNoF(indxModeNoF));
end
fracTracksDiffModeNoF = numDiffModeNoF / sum(numDiffModeNoF);

