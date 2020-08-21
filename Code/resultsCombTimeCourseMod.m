function [resSummaryInd] = resultsCombTimeCourseMod(dsSummary)
%RESULTSCOMBTIMECOURSE compiles the results of multiple movies from multiple timecourses into one timecourse
%
%SYNOPSIS [resSummaryComb,resSummaryInd] = resultsCombTimeCourseFunction(dsSummary,timeList,timeListComb,timeAbsOrRel)
%
%INPUT  dsSummary   : Structure array storing various results for each
%                     movie, ordered based on their time (corresponding
%                     to caseTimeList). See 'resultsIndTimeCourse'.
%       timeList    : 2-column vector indicating the movie times. Column 1
%                     shows absolute time, Column 2 shows relative time.
%       timeListComb: Column vector with timecourse time points for
%                     combined dataset.
%       timeAbsOrRel: 'abs' or 'rel' to use, respectively, absolute or
%                     relative time to align and group the different datasets.
%                     The absolute and relative times are stored in the 
%                     results file of each dataset.
%                     Optional. Default: 'abs'.
%    
%OUTPUT resSummaryComb: Structure with results summary for combined
%                       timecourse. Contains the fields: 
%           .numAbsClass: Absolute number of particles in the various
%                         motion classes.
%                         Row = time points in time course (see timeList).
%                         Columns = immobile, confined, free, directed,
%                         undetermined, determined, total.
%           .numNorm0Class: Normalized number of particles in the various
%                         motion classes, such that the mean of the first 3
%                         timepoints = 1.
%                         Rows and columns as numAbsClass.
%           .probClass  : Probability of the various motion classes.
%                         Rows as above.
%                         Columns = immobile, confined, free, directed (all
%                         relative to determined); and determined relative
%                         to total.
%           .diffCoefClass: Mean diffusion coefficient in the various
%                         motion classes.
%                         Rows as above.
%                         Columns = immobile, confined, free, directed,
%                         undetermined.
%           .confRadClass: Mean confinement radius in the various
%                         motion classes.
%                         Rows and columns as above.
%           .ampClass   : Mean amplitude of particles in the various motion
%                         classes.
%                         Rows and columns as above.
%           .ampNormClass:Mean normalized amplitude of particles in the
%                         various motion classes.
%                         Rows and columns as above.
%           .ampStatsF20: Amplitude statistics for particles in first 20
%                         frames.
%                         Rows as above.
%                         Columns = mean, first mode mean, first mode std,
%                         first mode fraction, number of modes, normalized
%                         mean.
%           .ampStatsL20: Amplitude statistics for particles in last 20
%                         frames.
%                         Rows and columns as above.
%           .rateMS     : Rate of merging and rate of splitting per
%                         feature (per time unit per number of features).
%                         Columns: merging, splitting.
%           .timeList   : List of time points in combined time course. Same
%                         as input timeListComb.
%           .timeAbsOrRel:'abs' or 'rel' to indicate, respectively, whether
%                         absolute or relative time was used to align and
%                         combine datasets.
%       resSummaryInd : Structure array equivalent to resSummaryComb but
%                       for individual timecourses.
%
%Khuloud Jaqaman, March 2015
%modified from resultsCombTimeCourse
%Tae H Kim, June 2015
%
%Further modifications: Khuloud Jaqaman, December 2017.
%Note: All diffusion mode analysis compilation and plotting makes the
%assumption that there is at most 10 modes. This is a hack, and should be
%fixed in the future.
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

%% Input

%get number of datasets
%numDS = 1;

%% Individual dataset results

%reserve memory
%{
[numAbsClassInd,numNorm0ClassInd,probClassInd,...
    diffCoefClassInd,confRadClassInd,ampClassInd,ampNormClassInd,...
    ampStatsF20Ind,ampStatsL20Ind,rateMSInd] = deal(cell(numDS,1));
%}
for iDS = 1
    
    %get number of time points in this movielist
    nTP = length(dsSummary);
    
    %tracks
    tracks = {dsSummary.tracks}';
        
    %determine if any time points are empty
    existTP = zeros(nTP,1);
    for iTP = 1 : nTP
        existTP(iTP) = ~isempty(tracks{iTP});
    end
    goodTP = find(existTP);
    numGoodTP = length(goodTP);
    missingAllTP = (numGoodTP == 0);
    
    %absolute numbers of molecules in various motion classes from MSS
    %analysis
    %columns are: 1 imm, 2 conf, 3 free, 4 dir, 5 undet, 6 det, 7 tot, 8 imm+conf
    %rows are for different timepoints, as listed in corresponding timeList
    numAbsClassInd = NaN(nTP,8);
    if ~missingAllTP
        diffSummary = vertcat(dsSummary.diffSummary); %#ok<NASGU>
        numAbsClassCurrent = [...
            catStruct(1,'diffSummary.probMotionType(6,3)') ...
            catStruct(1,'diffSummary.probMotionType(7,3)') ...
            catStruct(1,'diffSummary.probMotionType(8,3)') ...
            catStruct(1,'diffSummary.probMotionType(9,3)') ...
            catStruct(1,'diffSummary.probMotionType(11,3)')];
        numAbsClassInd(goodTP,:) = [numAbsClassCurrent sum(numAbsClassCurrent(:,1:4),2) ...
            sum(numAbsClassCurrent(:,1:5),2) sum(numAbsClassCurrent(:,1:2),2)];
    end
        
    %absolute numbers of molecules in various motion modes from diffusion
    %mode analysis
    %columns 1 to 10 are for modes 1 to 10.
    %column 11 is for tracks with undeteremined mode
    %column 12 is for all tracks with determined mode
    %column 13 is for total tracks
    numAbsModeInd = NaN(nTP,13);
    if ~missingAllTP
        diffModeSummary = vertcat(dsSummary.diffModeSummary); %#ok<NASGU>
        numAbsModeCurrent = catStruct(2,'diffModeSummary.probMotionMode(:,2)')';
        numModes = size(numAbsModeCurrent,2) - 1;
        if numModes < 10
            numAbsModeCurrent = [numAbsModeCurrent(:,1:numModes) NaN(numGoodTP,10-numModes) numAbsModeCurrent(:,numModes+1)];
        else
            warning('Number of modes exceeds 10; modes higher than 10 will be ignored')
            numAbsModeCurrent = numAbsModeCurrent(:,[1:10 end]);
        end
        numAbsModeInd(goodTP,:) = [numAbsModeCurrent nansum(numAbsModeCurrent(:,1:10),2) ...
            nansum(numAbsModeCurrent(:,1:11),2)];
    else
        numModes = 10;
    end
    
    %absolute probabilities of various motion classes from MSS analysis
    %columns are: 1 imm, 2 conf, 3 free, 4 dir (all relative to det), 
    %5 det (relative to tot), 6 imm+conf (relative to det)
    probAbsClassInd = [numAbsClassInd(:,1:4)./repmat(numAbsClassInd(:,6),1,4) ...
        numAbsClassInd(:,6)./numAbsClassInd(:,7) ...
        numAbsClassInd(:,8)./numAbsClassInd(:,6)];
    
    %absolute probabilities of various motion modes from diffusion mode
    %analysis
    %columns 1 to 10 are for modes 1 to 10 (relative to all determined)
    %column 11 is for all determined (relative to total)
    probAbsModeInd = [numAbsModeInd(:,1:10)./repmat(numAbsModeInd(:,12),1,10) ...
        numAbsModeInd(:,12)./numAbsModeInd(:,13)];
    
    %cell area in order to get densities
    cellAreaInd = NaN(nTP,1);
    if ~missingAllTP
        cellAreaInd(goodTP) = vertcat(dsSummary.cellArea);
    end
    
    %densities of molecules in various motion classes from MSS analysis
    %columns and rows same as absolute numbers
    densityAbsClassInd = numAbsClassInd ./ repmat(cellAreaInd,1,8);
    
    %densities of molecules in various motion modes from diffusion mode
    %analysis
    %columns and rows same as absolute numbers
    densityAbsModeInd = numAbsModeInd ./ repmat(cellAreaInd,1,13);
    
    %diffusion coefficient in various motion classes from MSS analysis
    %columns are: 1 imm, 2 conf, 3 free, 4 dir, 5 undet
    diffCoefClassInd = NaN(nTP,5);
    if ~missingAllTP
        diffCoefClassInd(goodTP,:) = (horzcat(dsSummary.diffCoefMeanPerClass))';
    end
    
    %confinement radius in various motion classes from MSS analysis
    %columns are: 1 imm, 2 conf, 3 free, 4 dir, 5 undet
    confRadClassInd = NaN(nTP,5);
    if ~missingAllTP
        confRadClassInd(goodTP,:) = (horzcat(dsSummary.confRadMeanPerClass))';
    end

    %diffusion coefficient in various motion modes from diffusion mode
    %analysis
    %columns are: modes 1 to N + undetermined
    diffCoefModeInd = NaN(nTP,numModes+1);
    if ~missingAllTP
        diffCoefModeInd(goodTP,:) = (horzcat(dsSummary.diffCoefMeanPerMode))';
    end
    
    %frame-to-frame mean square displacement in various motion modes from
    %diffusion mode analysis
    %columns are: modes 1 to N + undetermined
    f2fMeanSqDispModeInd = NaN(nTP,numModes+1);
    if ~missingAllTP
        f2fMeanSqDispModeInd(goodTP,:) = (horzcat(dsSummary.f2fMeanSqDispPerMode))';
    end
    
    %mean positional std in various motion modes from diffusion mode
    %analysis
    %columns are: modes 1 to N + undetermined
    meanPosStdModeInd = NaN(nTP,numModes+1);
    if ~missingAllTP
        meanPosStdModeInd(goodTP,:) = (horzcat(dsSummary.meanPosStdPerMode))';
    end
    
    %amplitude in various motion classes from MSS analysis
    %columns are: 1 imm, 2 conf, 3 free, 4 dir, 5 undet
    [ampClassInd,ampNormClassInd] = deal(NaN(nTP,5));
    if ~missingAllTP
        tmpCollect = (horzcat(dsSummary.ampMeanPerClass))';
        ampClassInd(goodTP,:) = tmpCollect(1:2:end,:);
        ampNormClassInd(goodTP,:) = tmpCollect(2:2:end,:);
    end
    
    %amplitude in various motion modes from diffusion mode analysis
    %columns are: modes 1 to N + undetermined
    [ampModeInd,ampNormModeInd] = deal(NaN(nTP,numModes+1));
    if ~missingAllTP
        tmpCollect = (horzcat(dsSummary.ampMeanPerMode))';
        ampModeInd(goodTP,:) = tmpCollect(1:2:end,:);
        ampNormModeInd(goodTP,:) = tmpCollect(2:2:end,:);
    end
    
    %patch up all variables to make 10 modes
    if numModes < 10
        diffCoefModeInd = [diffCoefModeInd(:,1:numModes) NaN(nTP,10-numModes) diffCoefModeInd(:,numModes+1)];
        f2fMeanSqDispModeInd = [f2fMeanSqDispModeInd(:,1:numModes) NaN(nTP,10-numModes) f2fMeanSqDispModeInd(:,numModes+1)];
        meanPosStdModeInd = [meanPosStdModeInd(:,1:numModes) NaN(nTP,10-numModes) meanPosStdModeInd(:,numModes+1)];
        ampModeInd = [ampModeInd(:,1:numModes) NaN(nTP,10-numModes) ampModeInd(:,numModes+1)];
        ampNormModeInd = [ampNormModeInd(:,1:numModes) NaN(nTP,10-numModes) ampNormModeInd(:,numModes+1)];
   else
        diffCoefModeInd = diffCoefModeInd(:,[1:10 end]);
        f2fMeanSqDispModeInd = f2fMeanSqDispModeInd(:,[1:10 end]);
        meanPosStdModeInd = meanPosStdModeInd(:,[1:10 end]);
        ampModeInd = ampModeInd(:,[1:10 end]);
        ampNormModeInd = ampNormModeInd(:,[1:10 end]);
    end
    
    %diffusion mode decomposition
    [numDiffModeInd,numDiffModeControlInd] = deal(NaN(nTP,1));
    [modeDiffCoefInd,modeFractionInd] = deal(NaN(nTP,10));
    if ~missingAllTP
        numDiffModeInd(goodTP) = catStruct(1,'dsSummary.diffModeDecomposition.numMode');
        numDiffModeControlInd(goodTP) = catStruct(1,'dsSummary.diffModeDecomposition.numModeControl');
        for iMode = 1 : 10
            indxGood = find(numDiffModeInd>=iMode);
            dsSummaryGood = dsSummary(indxGood); %#ok<NASGU>
            modeDiffCoefInd(indxGood,iMode) = catStruct(1,['dsSummaryGood.diffModeDecomposition.modeParam(' num2str(iMode) ',1)']);
            modeFractionInd(indxGood,iMode) = catStruct(1,['dsSummaryGood.diffModeDecomposition.modeParam(' num2str(iMode) ',2)']);
        end
    end
    
    %amplitude statistics in first 20 frames
    %columns are: 1: mean, 2-4: modes 1-3 mean, 5-7: modes 1-3 std, 8-10:
    %modes 1 2 >2 fractions, 11: number of modes, 12: normalized mean
    ampStatsF20Ind = NaN(nTP,12);
    if ~missingAllTP
        ampStatsF20Ind(goodTP,:) = vertcat(dsSummary.ampStatsF20);
    end
    
    %amplitude statistics in last 20 frames
    %columns are same as above
    ampStatsL20Ind = NaN(nTP,12);
    if ~missingAllTP
        ampStatsL20Ind(goodTP,:) = vertcat(dsSummary.ampStatsL20);
    end
    
    %amplitude statistics in first frame from detection
    %columns are same as above
    ampStatsF01Ind = NaN(nTP,12);
    if ~missingAllTP
        ampStatsF01Ind(goodTP,:) = vertcat(dsSummary.ampStatsF01);
    end
    
    %rate of merging and rate of splitting
    rateMSInd = NaN(nTP,2);
    if ~missingAllTP
        tmp = vertcat(dsSummary.statsMS);
        rateMSInd(goodTP,:) = tmp(:,8:9);
    end
    
    %merging and splitting time information
    msTimeInfoInd = NaN(nTP,5);
    if ~missingAllTP
        msTimeInfoInd(goodTP,:) = vertcat(dsSummary.msTimeInfo);
    end
    
end

%% Combine dataset results
%{
%do direct combination for all variables except for
%normalized number of particles in various classes
[numAbsClassComb.msn,numAbsClassComb.sample] = combineWithTimeAlign(numAbsClassInd,timeList,timeListComb,timeIndxComb);
[probClassComb.msn,probClassComb.sample] = combineWithTimeAlign(probClassInd,timeList,timeListComb,timeIndxComb);
[diffCoefClassComb.msn,diffCoefClassComb.sample] = combineWithTimeAlign(diffCoefClassInd,timeList,timeListComb,timeIndxComb);
[confRadClassComb.msn,confRadClassComb.sample] = combineWithTimeAlign(confRadClassInd,timeList,timeListComb,timeIndxComb);
[ampClassComb.msn,ampClassComb.sample] = combineWithTimeAlign(ampClassInd,timeList,timeListComb,timeIndxComb);
[ampNormClassComb.msn,ampNormClassComb.sample] = combineWithTimeAlign(ampNormClassInd,timeList,timeListComb,timeIndxComb);
[ampStatsF20Comb.msn,ampStatsF20Comb.sample] = combineWithTimeAlign(ampStatsF20Ind,timeList,timeListComb,timeIndxComb);
[ampStatsL20Comb.msn,ampStatsL20Comb.sample] = combineWithTimeAlign(ampStatsL20Ind,timeList,timeListComb,timeIndxComb);
[rateMSComb.msn,rateMSComb.sample] = combineWithTimeAlign(rateMSInd,timeList,timeListComb,timeIndxComb);

%for normalized number of particles, normalize by mean value at time 0 for
%combined datasets
tmp = numAbsClassComb.sample;
tmp = tmp ./ repmat(numAbsClassComb.msn(1,:,1),[size(tmp,1) 1 size(tmp,3)]);
numNorm0ClassComb.sample = tmp;
numNorm0ClassComb.msn = numAbsClassComb.msn;
numNorm0ClassComb.msn(:,:,1) = nanmean(tmp,3);
numNorm0ClassComb.msn(:,:,2) = nanstd(tmp,[],3);
%}
%% Output
%{
%Replace all nan with 0
numAbsClassInd(isnan(numAbsClassInd)) = 0;
numNorm0ClassInd(isnan(numNorm0ClassInd)) = 0;
probClassInd(isnan(probClassInd)) = 0;
diffCoefClassInd(isnan(diffCoefClassInd)) = 0;
confRadClassInd(isnan(confRadClassInd)) = 0;
ampClassInd(isnan(ampClassInd)) = 0;
ampNormClassInd(isnan(ampNormClassInd)) = 0;
ampStatsF20Ind(isnan(ampStatsF20Ind)) = 0;
ampStatsL20Ind(isnan(ampStatsL20Ind)) = 0;
%}

% resSummaryInd = struct('numAbsClass',numAbsClassInd,'numNorm0Class',numNorm0ClassInd,...
%     'densityAbsClass',densityAbsClassInd,'densityNorm0Class',densityNorm0ClassInd,...
%     'probAbsClass',probAbsClassInd,'probNorm0Class',probNorm0ClassInd,...
%     'diffCoefClass',diffCoefClassInd,'confRadClass',confRadClassInd,...
%     'ampClass',ampClassInd,'ampNormClass',ampNormClassInd,...
%     'ampStatsF20',ampStatsF20Ind,'ampStatsL20',ampStatsL20Ind,'ampStatsF01',ampStatsF01Ind,...
%     'rateMS',rateMSInd,'msTimeInfo',msTimeInfoInd);

resSummaryInd = struct('numAbsClass',numAbsClassInd,'numAbsMode',numAbsModeInd,...
    'densityAbsClass',densityAbsClassInd,'densityAbsMode',densityAbsModeInd,...
    'probAbsClass',probAbsClassInd,'probAbsMode',probAbsModeInd,...
    'diffCoefClass',diffCoefClassInd,'confRadClass',confRadClassInd,...
    'diffCoefMode',diffCoefModeInd,'f2fMeanSqDispMode',f2fMeanSqDispModeInd,...
    'meanPosStdMode',meanPosStdModeInd,...
    'ampClass',ampClassInd,'ampNormClass',ampNormClassInd,...
    'ampMode',ampModeInd,'ampNormMode',ampNormModeInd,...
    'ampStatsF20',ampStatsF20Ind,'ampStatsL20',ampStatsL20Ind,'ampStatsF01',ampStatsF01Ind,...
    'rateMS',rateMSInd,'msTimeInfo',msTimeInfoInd,...
    'numDiffMode',[numDiffModeInd numDiffModeControlInd],...
    'modeDiffCoef',modeDiffCoefInd,'modeFraction',modeFractionInd); 
resSummaryInd.tracks = tracks;

%% ~~~ the end ~~~
