function [ resSummary ] = resultsIndTimeCoursePerMovie( MD, saveFile, channels )
%resultsIndTimeCoursePerMovie Subroutine of resultsIndTimeCourse
%
% INPUT
% MD        : MovieData object or string for MovieData .mat file
% saveFile  : File to which to save resSummary (optional, default: no save)
% channels : (optional) What channels to use, default:
%                       1:length(MD.channels_)
% See also resultsIndTimeCourse
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

% Based on resultsIndTimeCourse by Khuloud Jaqaman, March 2015
% Mark Kittisopikul, October 2015

if(ischar(MD))
    MD = MovieData.load(MD);
end
if(nargin < 3 || isempty(channels))
    channels = 1 : length(MD.channels_);
    curChannels = channels;
else
    curChannels = channels(:)';
    curChannels = curChannels(curChannels <= length(MD.channels_));
end

% channels could be empty or exceed the number of channels in this movie
resSummary(1,max(channels)) = timeCourseAnalysis.util.emptyResSummaryStruct;

progressTextMultiple('channel', length(curChannels));

for iC = curChannels
    
    iProcTracks = MD.getProcessIndex('TrackingProcess',1,0); %tracking process
    if isempty(iProcTracks)
        error([MD.movieDataPath_ ' : Tracking Process Missing']);
    end
    trackData = load(MD.processes_{iProcTracks}.outFilePaths_{iC});
    tracks0 = trackData.tracksFinal;
    
    %continue only if there are tracks in this movie, otherwise return empty
    if ~isempty(tracks0)
        
        iProcDiff = MD.getProcessIndex('MotionAnalysisProcess',1,0); %diffusion analysis and tracks used for diffusion analysis
        if isempty(iProcDiff)
            error([MD.movieDataPath_ ' : Motion Analysis Process Missing']);
        end
        motionAnalysis = load(MD.processes_{iProcDiff}.outFilePaths_{iC});
        
        iProcMask = MD.getProcessIndex('ImportCellMaskProcess',1,0); %cell mask
        if ~isempty(iProcMask)
            % TODO: Is there a different mask for the channel?
            mask = imread(fullfile(MD.processes_{iProcMask}.funParams_.OutputDirectory,'cellMask_channel_1.tif'));
            cellArea = sum(mask(:));
        else
            mask = [];
            cellArea = 512^2; %hard-code for now - look for better solution given the peculiarities of our simultaneous 2-color imaging
        end
        
        iProcDet = MD.getProcessIndex('DetectionProcess',1,0); %detection - if code makes it to here, then detection has been done
        detectionOutput = load(MD.processes_{iProcDet}.outFilePaths_{iC});
        movieInfoFramesAll = detectionOutput.movieInfo;
        
        %get number of frames and first frame to use for analysis
        numFrames = length(movieInfoFramesAll);
        numMol = zeros(numFrames,1);
        for iFrame = 1 : numFrames
            numMol(iFrame) = size(movieInfoFramesAll(iFrame).xCoord,1);
        end
        firstFrameWithMol = find(numMol > 0, 1, 'first');
        if firstFrameWithMol > 1
            firstFrameToUse = min(firstFrameWithMol+3,numFrames); %this is a hack to handle Aparajita's data - be careful, it might cause problems for other scenarios
        else
            firstFrameToUse = firstFrameWithMol;
        end
        lastFrameToUse = min(firstFrameToUse+4,numFrames);
        
        iProcDiffModes = MD.getProcessIndex('DiffusionModeAnalysisProcess',1,0); %diffusion mode analysis
        if isempty(iProcDiffModes)
            error([MD.movieDataPath_ ' : Diffusion Mode Analysis Process Missing']);
        end
        diffModeAnalysis = load(MD.processes_{iProcDiffModes}.outFilePaths_{iC});
        try
            numModes = size(diffModeAnalysis.diffModeAnalysisRes(1).summary.probMotionMode,1) - 1; %number of modes, needed in various instances later
        catch
            numModes = NaN;
        end
        
        %limit analysis to tracks in mask if supplied
        if ~isempty(mask) && any(mask(:)==0)
            
            %do this for overall tracks
            numTracks = length(tracks0);
            keepTrack = ones(numTracks,1);
            for iTrack = 1 : numTracks
                xCoord = tracks0(iTrack).tracksCoordAmpCG(:,1:8:end);
                yCoord = tracks0(iTrack).tracksCoordAmpCG(:,2:8:end);
                meanPosX = round(nanmean(xCoord(:)));
                meanPosY = round(nanmean(yCoord(:)));
                keepTrack(iTrack) = mask(meanPosY,meanPosX);
            end
            tracks0 = tracks0(find(keepTrack));
            
        end
        
        %continue only if there are tracks in mask, otherwise return empty
        if ~isempty(tracks0)
            
            if ~isempty(mask) && any(mask(:)==0)
                
                %do this for tracks used in conventional diffusion analysis
                numTracks = length(motionAnalysis.tracks);
                keepTrack = ones(numTracks,1);
                for iTrack = 1 : numTracks
                    xCoord = motionAnalysis.tracks(iTrack).tracksCoordAmpCG(:,1:8:end);
                    yCoord = motionAnalysis.tracks(iTrack).tracksCoordAmpCG(:,2:8:end);
                    meanPosX = round(nanmean(xCoord(:)));
                    meanPosY = round(nanmean(yCoord(:)));
                    keepTrack(iTrack) = mask(meanPosY,meanPosX);
                end
                indxKeep = find(keepTrack);
                motionAnalysis.tracks = motionAnalysis.tracks(indxKeep);
                motionAnalysis.diffAnalysisRes = motionAnalysis.diffAnalysisRes(indxKeep);
                
                %redo conventional diffusion analysis summary
                if length(indxKeep) < numTracks
                    minTrackLen = 5;
                    probDim = 2;
                    extractType = 1;
                    [probMotionType,motionChar] = summarizeDiffAnRes(...
                        motionAnalysis.tracks,minTrackLen,probDim,...
                        motionAnalysis.diffAnalysisRes,extractType);
                    motionAnalysis.diffAnalysisRes(1).summary.probMotionType = probMotionType;
                    motionAnalysis.diffAnalysisRes(1).summary.motionChar = motionChar;
                end
                
                %do this for tracks used in diffusion mode analysis
                numTracks = length(diffModeAnalysis.tracks);
                keepTrack = ones(numTracks,1);
                for iTrack = 1 : numTracks
                    xCoord = diffModeAnalysis.tracks(iTrack).tracksCoordAmpCG(:,1:8:end);
                    yCoord = diffModeAnalysis.tracks(iTrack).tracksCoordAmpCG(:,2:8:end);
                    meanPosX = round(nanmean(xCoord(:)));
                    meanPosY = round(nanmean(yCoord(:)));
                    keepTrack(iTrack) = mask(meanPosY,meanPosX);
                end
                indxKeep = find(keepTrack);
                diffModeAnalysis.tracks = diffModeAnalysis.tracks(indxKeep);
                diffModeAnalysis.diffModeAnalysisRes = diffModeAnalysis.diffModeAnalysisRes(indxKeep);
                
                %redo diffusion mode analysis summary
                if length(indxKeep) < numTracks
                    minTrackLen = 5;
                    probDim = 2;
                    extractType = 1;
                    [probMotionMode,modeMotionChar] = summarizeDiffModeRes(...
                        diffModeAnalysis.tracks,diffModeAnalysis.diffModeAnalysisRes,...
                        numModes,minTrackLen,probDim,extractType);
                    diffModeAnalysis.diffModeAnalysisRes(1).summary.probMotionMode = probMotionMode;
                    diffModeAnalysis.diffModeAnalysisRes(1).summary.modeMotionChar = modeMotionChar;
                end
                
            end
            
            %conventional diffusion analysis summary
            diffSummary = motionAnalysis.diffAnalysisRes(1).summary;
            
            %conventional diffusion analysis per track
            trajClass = vertcat(motionAnalysis.diffAnalysisRes.classification);
            trajClass = trajClass(:,2);
            trajDiffCoef = catStruct(1,'motionAnalysis.diffAnalysisRes.fullDim.normDiffCoef');
            trajConfRad = catStruct(1,'motionAnalysis.diffAnalysisRes.confRadInfo.confRadius');
            
            %diffusion mode decomposition
            %nothing to do here, just copy results directly into resSummary
            
            %diffusion mode analysis summary
            diffModeSummary = diffModeAnalysis.diffModeAnalysisRes(1).summary;
            
            %diffusion mode analysis per track
            trajMode = vertcat(diffModeAnalysis.diffModeAnalysisRes.diffMode);
            trajDiffCoefMode = vertcat(diffModeAnalysis.diffModeAnalysisRes.diffCoef);
            trajMSDF2F = vertcat(diffModeAnalysis.diffModeAnalysisRes.msdF2F);
            trajMeanPosStd = vertcat(diffModeAnalysis.diffModeAnalysisRes.meanPosStd);
            %     trajDiffRadius = vertcat(diffModeAnalysis.diffModeAnalysisRes.diffRadius);
            %     trajLifetime = vertcat(diffModeAnalysis.diffModeAnalysisRes.lifetime);
            
            %amplitude matrix
            tracksMat = convStruct2MatIgnoreMS(tracks0);
            trackSEL = getTrackSEL(tracksMat);
            indxKeepLater = find(trackSEL(:,3)>=3);
            ampMatTmp = tracksMat(:,4:8:end);
            [ampNumRow,ampNumCol] = size(ampMatTmp);
            ampMat = NaN(ampNumRow,numFrames);
            ampMat(:,1:ampNumCol) = ampMatTmp;
            
            %amplitude per track
            %KJ, 171219: get amplitudes only in first 5 frames of movie
            %(instead of all throughout), to minimize effect of photobleaching.
            %This means not all tracks will get an amplitude
            ampMeanPerTraj = nanmean(ampMat(:,firstFrameToUse:lastFrameToUse),2);
            
            %average properties per motion class
            [ampMeanPerClass,diffCoefMeanPerClass,confRadMeanPerClass] = deal(NaN(5,1));
            for i = 1 : 2
                indxClass = find(trajClass == (i-1));
                ampMeanPerClass(i) = nanmean(ampMeanPerTraj(indxClass));
                diffCoefMeanPerClass(i) = mean(trajDiffCoef(indxClass));
                confRadMeanPerClass(i) = mean(trajConfRad(indxClass));
            end
            for i = 3 : 4
                indxClass = find(trajClass == (i-1));
                ampMeanPerClass(i) = nanmean(ampMeanPerTraj(indxClass));
                diffCoefMeanPerClass(i) = mean(trajDiffCoef(indxClass));
            end
            ampMeanPerClass(5) = nanmean(ampMeanPerTraj(isnan(trajClass)));
            
            %average properties per diffusion mode
            [ampMeanPerMode,diffCoefMeanPerMode,f2fMeanSqDispPerMode,meanPosStdPerMode] = deal(NaN(numModes+1,1));
            for i = 1 : numModes
                indxMode = find(trajMode == i);
                ampMeanPerMode(i) = nanmean(ampMeanPerTraj(indxMode));
                diffCoefMeanPerMode(i) = mean(trajDiffCoefMode(indxMode));
                f2fMeanSqDispPerMode(i) = mean(trajMSDF2F(indxMode));
                meanPosStdPerMode(i) = mean(trajMeanPosStd(indxMode));
                %         diffRadiusMeanPerMode(i) = mean(trajDiffRadius(indxMode));
                %         lifePerMode(i) = mean(trajLifetime(indxMode));
            end
            %unclassified tracks
            indxMode = find(isnan(trajMode));
            ampMeanPerMode(end) = nanmean(ampMeanPerTraj(indxMode));
            diffCoefMeanPerMode(end) = mean(trajDiffCoefMode(indxMode));
            f2fMeanSqDispPerMode(end) = mean(trajMSDF2F(indxMode));
            meanPosStdPerMode(end) = mean(trajMeanPosStd(indxMode));
            
            %amplitude statistics in first 5 frames
%             indxMode1 = find(trajMode == 1);
%             ampVec = ampMat(intersect(indxKeepLater,indxMode1),firstFrameToUse:lastFrameToUse);
            ampVec = ampMat(indxKeepLater,firstFrameToUse:lastFrameToUse);
            ampVec = ampVec(~isnan(ampVec));
            maxGaussians = min(4,numel(ampVec)-3);
            try
                if maxGaussians > 0
                    %                 [~,~,modeParam] = fitHistWithGaussians(ampVec,1,0,3,0,[1 maxGaussians],2,[],1,[],0);
                    [~,~,modeParam] = fitHistWithGaussians(ampVec,1,-2,3,0,[1 maxGaussians],2,[],1,[],0.1);
                    numMode = size(modeParam,1);
                    modeParamMean = modeParam(:,1)';
                    modeParamStd  = modeParam(:,2)';
                    ampModeMean = exp( modeParamMean + modeParamStd.^2/2 );
                    ampModeMean = [ampModeMean(1:min(3,numMode)) NaN(1,max(0,3-numMode))];
                    ampModeStd  = sqrt( exp( modeParamStd.^2 + 2*modeParamMean ) .* ( exp( modeParamStd.^2 )-1 ) );
                    ampModeStd  = [ampModeStd(1:min(3,numMode)) NaN(1,max(0,3-numMode))];
                    ampModeFrac = modeParam(:,4)'/length(ampVec);
                    ampModeFrac = [ampModeFrac zeros(1,abs(2-numMode))]; %#ok<AGROW>
                    ampModeFrac = [ampModeFrac(1:2) 1-sum(ampModeFrac(1:2))];
                    ampStatsF20 = [mean(ampVec) ampModeMean ampModeStd ampModeFrac numMode];
                else
                    ampStatsF20 = [mean(ampVec) NaN(1,9) 0];
                end
            catch
                ampStatsF20 = [mean(ampVec) NaN(1,9) 0];
            end
            
            %amplitude statistics in last 5 frames
            ampVec = ampMat(indxKeepLater,max(1,end-4):end);
            ampVec = ampVec(~isnan(ampVec));
            maxGaussians = min(4,numel(ampVec)-3);
            try
                if maxGaussians > 0
                    %                 [~,~,modeParam] = fitHistWithGaussians(ampVec,1,0,3,0,[1 maxGaussians],2,[],1,[],0);
                    [~,~,modeParam] = fitHistWithGaussians(ampVec,1,-2,3,0,[1 maxGaussians],2,[],1,[],0.1);
                    numMode = size(modeParam,1);
                    modeParamMean = modeParam(:,1)';
                    modeParamStd  = modeParam(:,2)';
                    ampModeMean = exp( modeParamMean + modeParamStd.^2/2 );
                    ampModeMean = [ampModeMean(1:min(3,numMode)) NaN(1,max(0,3-numMode))];
                    ampModeStd  = sqrt( exp( modeParamStd.^2 + 2*modeParamMean ) .* ( exp( modeParamStd.^2 )-1 ) );
                    ampModeStd  = [ampModeStd(1:min(3,numMode)) NaN(1,max(0,3-numMode))];
                    ampModeFrac = modeParam(:,4)'/length(ampVec);
                    ampModeFrac = [ampModeFrac zeros(1,abs(2-numMode))]; %#ok<AGROW>
                    ampModeFrac = [ampModeFrac(1:2) 1-sum(ampModeFrac(1:2))];
                    ampStatsL20 = [mean(ampVec) ampModeMean ampModeStd ampModeFrac numMode];
                else
                    ampStatsL20 = [mean(ampVec) NaN(1,9) 0];
                end
            catch
                ampStatsL20 = [mean(ampVec) NaN(1,9) 0];
            end
            
            %amplitude statistics in first frame directly from detection
            ampVec = movieInfoFramesAll(firstFrameToUse).amp;
            if ~isempty(ampVec)
                ampVec = ampVec(:,1);
                maxGaussians = min(4,numel(ampVec)-3);
            else
                maxGaussians = 0;
            end
            try
                if maxGaussians > 0
                    %                 [~,~,modeParam] = fitHistWithGaussians(ampVec,1,0,3,0,[1 maxGaussians],2,[],1,[],0);
                    [~,~,modeParam] = fitHistWithGaussians(ampVec,1,-2,3,0,[1 maxGaussians],2,[],1,[],0.1);
                    numMode = size(modeParam,1);
                    modeParamMean = modeParam(:,1)';
                    modeParamStd  = modeParam(:,2)';
                    ampModeMean = exp( modeParamMean + modeParamStd.^2/2 );
                    ampModeMean = [ampModeMean(1:min(3,numMode)) NaN(1,max(0,3-numMode))];
                    ampModeStd  = sqrt( exp( modeParamStd.^2 + 2*modeParamMean ) .* ( exp( modeParamStd.^2 )-1 ) );
                    ampModeStd  = [ampModeStd(1:min(3,numMode)) NaN(1,max(0,3-numMode))];
                    ampModeFrac = modeParam(:,4)'/length(ampVec);
                    ampModeFrac = [ampModeFrac zeros(1,abs(2-numMode))]; %#ok<AGROW>
                    ampModeFrac = [ampModeFrac(1:2) 1-sum(ampModeFrac(1:2))];
                    ampStatsF01 = [mean(ampVec) ampModeMean ampModeStd ampModeFrac numMode];
                else
                    ampStatsF01 = [mean(ampVec) NaN(1,9) 0];
                end
            catch
                ampStatsF01 = [mean(ampVec) NaN(1,9) 0];
            end
            
            %normalize mean amplitudes by first mode mean
            %for amp per class, amp per mode, and amp in first 5 frames, use mode of first 5 frames
            ampMeanPerClass = [ampMeanPerClass ampMeanPerClass/(ampStatsF20(3)/2)]; %#ok<AGROW>
            ampMeanPerMode = [ampMeanPerMode ampMeanPerMode/(ampStatsF20(3)/2)]; %#ok<AGROW>
            ampStatsF20 = [ampStatsF20 ampStatsF20(1)/(ampStatsF20(3)/2)]; %#ok<AGROW>
            %for last 5 frames, use mode of last 5 frames
            ampStatsL20 = [ampStatsL20 ampStatsL20(1)/(ampStatsL20(3)/2)]; %#ok<AGROW>
            %for first frame from detection directly, use mode of first frame
            ampStatsF01 = [ampStatsF01 ampStatsF01(1)/(ampStatsF01(3)/2)]; %#ok<AGROW>
            
            %merge and split statistics
            statsMS = calcStatsMS_noMotionInfo(tracks0,5,1,1);
            tmp = calcMergeSplitTimes(tracks0,5,[],1);
            tmp = tmp.all;
            msTimeInfo = [mean(tmp.timeMerge2Split) mean(tmp.timeSplit2MergeSelf) mean(tmp.timeSplit2MergeOther) mean(tmp.timeMerge2End) mean(tmp.timeStart2Split)];
            %         msTimeInfo = tmp.timeMerge2Split;
            
            %results for output
            %conventional diffusion analysis
            resSummary(1,iC).diffSummary = diffSummary;
            resSummary(1,iC).diffCoefMeanPerClass = diffCoefMeanPerClass;
            resSummary(1,iC).confRadMeanPerClass = confRadMeanPerClass;
            %diffusion mode analysis
            resSummary(1,iC).diffModeDecomposition = diffModeAnalysis.diffModeDecomposition;
            resSummary(1,iC).diffModeSummary = diffModeSummary;
            resSummary(1,iC).diffCoefMeanPerMode = diffCoefMeanPerMode;
            resSummary(1,iC).f2fMeanSqDispPerMode = f2fMeanSqDispPerMode;
            resSummary(1,iC).meanPosStdPerMode = meanPosStdPerMode;
            %amplitude per motion class or mode
            resSummary(1,iC).ampMeanPerClass = ampMeanPerClass;
            resSummary(1,iC).ampMeanPerMode = ampMeanPerMode;
            %overall amplitude decomposition
            resSummary(1,iC).ampStatsF20 = ampStatsF20;
            resSummary(1,iC).ampStatsL20 = ampStatsL20;
            resSummary(1,iC).ampStatsF01 = ampStatsF01;
            %merging and splitting
            resSummary(1,iC).statsMS = statsMS;
            resSummary(1,iC).msTimeInfo = msTimeInfo;
            %other
            resSummary(1,iC).cellArea = cellArea;
            resSummary(1,iC).tracks = tracks0;
            %     resSummary(1,iC).diffModeAnalysis = struct('modeParam',ones(1,4),'numMode',1,'modeParamControl',ones(1,4),'numModeControl',1);
            
        end %(if ~isempty(tracks0))
        
    end %(if ~isempty(tracks0))
    
    %Progress Counter
    progressTextMultiple();
    
end %(for iC = curChannels)

if(nargin > 1 && ischar(saveFile))
    try
        save(saveFile,'resSummary');
    catch err
        warning(['resultIndTimeCoursePerMovie: Could not save to file, ' saveFile]);
        disp(getReport(err));
    end
end

progressTextMultiple();

end

