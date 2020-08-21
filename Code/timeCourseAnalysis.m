function [] = timeCourseAnalysis(CMLs, outputDir, varargin)
%TimeCourseAnalysis of CombinedMovieList objects.
%In TimeCourseAnalysis, MLs in each CML are considered to be in similar
%condition and grouped together for plotting.
%
%SYNOPSIS function [] = timeCourseAnalysis(CMLs, outputDir, varargin)
%
%INPUT
%   CMLs        : (array of CombinedMovieList) objects or (cell of string)
%                 containingfull paths to CMLs
%   ouputDir    : (string) Directory where figures and analysis result will
%                 be stored
%   varargin    : name_value pairs for analysis parameter
%       'smoothingPara'         : Parameter used for smoothing spline fit.
%       'channels'              : Channel of MD to be analyzed. Default: 1.
%       'doNewAnalysis'         : (logical) True: always do new analysis even if
%                                 the analysis has already been done.
%                                 False: avoid doing the analysis again if
%                                 analysis has already been done and the
%                                 analysis parameters are identical.
%                                 Default: true.
%       'doPartitionAnalysis'   : (logical) To analyze trackPartitioning
%                                 process or not. Default: false.
%       'shiftPlotPositive'     : (logical) Determines whether or not whole
%                                 plots are shifted so that no negative
%                                 time values are present. In other words,
%                                 minimum time point is taken to be the
%                                 zero. Default: false. <don't use>.
%       'start2zero'            : (logical) Sets the average start time to
%                                 be zero, even if the align event is not
%                                 'start'.
%       'channelNames'          : (cellstr) Cell array of channel names.
%       'detectOutliers_k_sigma': (numeric, scalar) Dee detectOutliers.
%                                 Default: 4.
%       'aveInterval'           : (numeric, scalar) Interval for averaging
%                                 single cell behavior over time course.
%                                 Default: 3 min.
%
%
%Tae H Kim, July 2015
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

%convert CMLs if it's cell of strings
nCML = numel(CMLs);
if iscellstr(CMLs)
    directory_CML = CMLs;
    clear CMLs
    for iCML = 1:nCML
        fprintf('Loading CombinedMovieList %g/%g\n', iCML, nCML);
        CMLs(iCML) = CombinedMovieList.load(directory_CML{iCML});
    end
elseif isa(CMLs, 'CombinedMovieList')
    directory_CML = arrayfun(@(x) x.getFullPath(), CMLs, 'UniformOutput', false); %#ok<NASGU>
else
    error('CMLs must be class object CombinedMovieList');
end
%input parser
ip = inputParser;
ip.CaseSensitive = false;
% Extra (Unmatched) parameters forwarded to timeCourseAnalysis_StandAlone
ip.KeepUnmatched = true;
ip.StructExpand = true;
ip.addRequired('outputDir', @ischar);
ip.addParameter('channels', [], @isnumeric);
ip.addParameter('doPartition', false, @(x) isnumeric(x) || islogical(x));
ip.addParameter('doNewAnalysis', true, @(x) isnumeric(x) || islogical(x));
ip.addParameter('start2zero', false, @(x) islogical(x)||isnumeric(x));
ip.addParameter('channelNames', false, @(x) iscellstr(x));
% ip.addParameter('detectOutliers_k_sigma', 4, @(x) isnumeric(x));
% ip.addParameter('aveInterval', 3, @(x) isnumeric(x));
ip.parse(outputDir, varargin{:});
%for easier access to ip.Result variables
analysisPara = ip.Results;
%checks if CMLs are loaded or not
%AND determines how many ML there are total
nMLTot = 0;
for iCML = 1:nCML
    progressDisp = fprintf('Loading Combined Movie List %g/%g\n', iCML, nCML); %loading progress display
    % if not loads it
    if numel(CMLs(iCML).movieLists_) ~= numel(CMLs(iCML).movieListDirectory_)
        CMLs(iCML).sanityCheck();
    end
    %checks if all MD has necessary processes
    arrayfun(@(x) timeCourseAnalysis.MLCheck(x, analysisPara), CMLs(iCML).movieLists_);
    fprintf(repmat('\b', 1, progressDisp)); %loading progress display
    nMLTot = nMLTot + numel(CMLs(iCML).movieLists_);
end

%% Main Time Course Analysis (CML-level)

% Analyze all MovieData in parallel first
timeCourseAnalysis.analyzeMDsInParallel(CMLs,analysisPara.doNewAnalysis,analysisPara.channels);
%Progress Display
progressTextMultiple('Time Course Analysis', nMLTot);
%Using resultsIndTimeCourseMod.m to do basic analysis
%and extract time data and align
[summary, time, extra, startTime] = arrayfun(@(CML) timeCourseAnalysis.CMLAnalyze(CML,analysisPara), CMLs, 'UniformOutput', false);
startTime = [startTime{:}];

%% Format Change

% Break summary apart by channel
summary = cellfun(@(s) num2cell(s,1),summary,'UniformOutput',false);
% Should a nCML x nChannel cell array
% Made an assumption that all that CMLs have the same channel length here
summary = vertcat(summary{:});
%Using resultsCombTimeCourseMod.m to store data in correct format
summaryTmp = cell(size(summary));
summaryTmp(:,analysisPara.channels) = cellfun(@resultsCombTimeCourseMod, summary(:,analysisPara.channels), 'UniformOutput', false);
summary = summaryTmp;
% summary = cellfun(@resultsCombTimeCourseMod, summary, 'UniformOutput', false);

%adds time and name (because that's not in resultsCombTimeCourseMod.m)
%KJ: also specify which movieList each movie comes from. We are making the
%assumption that movies are in correct order within each movieList
%(from first to last)
if(isempty(analysisPara.channels))
    analysisPara.channels = 1:size(summary,2);
end
for iCML = 1:nCML
    nML = length(CMLs(iCML).movieLists_);
    mlIndex = [];
    for iML = 1 : nML
        nMD = length(CMLs(iCML).movieLists_(iML).movies_);
        mlIndex = [mlIndex; iML*ones(nMD,1)]; %#ok<AGROW>
    end
    for iChannel = analysisPara.channels
        summary{iCML,iChannel}.time = time{iCML};
        summary{iCML,iChannel}.name = analysisPara.channelNames{iChannel};
        if(~isempty(CMLs(iCML).name_))
            if(~isempty(summary{iCML,iChannel}.name))
                summary{iCML,iChannel}.name = [ '/' summary{iCML,iChannel}.name];
            end
            summary{iCML,iChannel}.name = [ CMLs(iCML).name_ summary{iCML,iChannel}.name];
        end
        summary{iCML,iChannel}.mlIndex = mlIndex;
    end
end

%% Calculate normalized molecular properties (normalized by first 3 movies in each ML)

%after sorting, calculate some normalized properties
summaryTmp = cell(size(summary));
summaryTmp(:,analysisPara.channels) = cellfun(@resultsNormTimeCourse, summary(:,analysisPara.channels), 'UniformOutput', false);
summary = summaryTmp;
% summary = cellfun(@resultsNormTimeCourse, summary, 'UniformOutput', false);

%% Linearize indexing (HACK, FIX)

% essentially let the rest of the analysis use linear indexing to access
% summary
nCML = numel(summary);

%% Partitioning analysis, which is no longer used, and is thus most likely outdated

%adds extra analysis if applicable
if analysisPara.doPartition
    for iCML = 1:nCML
        summary{iCML}.chemEnergy = vertcat(extra{iCML}.chemEnergy);
        summary{iCML}.locFreq = vertcat(extra{iCML}.locFreq);
        summary{iCML}.delocFreq = vertcat(extra{iCML}.delocFreq);
        summary{iCML}.eqCond = vertcat(extra{iCML}.eqCond);
    end
end

%% Sort

%order data from earliest time point to latest
goodCML = zeros(nCML,1);
for iCML = 1:nCML
    if ~isempty(summary{iCML})
        [~, sortIndx] = sort(summary{iCML}.time);
        % Name should not be sorted
        name = summary{iCML}.name;
        summary{iCML} = structfun(@(x) x(sortIndx,:),summary{iCML},'UniformOutput',false,'ErrorHandler',@(~,x) x);
        summary{iCML}.name = name;
        goodCML(iCML) = 1;
    end
end
goodCML = find(goodCML);

%% Shift time

%% ADD MULTIPLE CHANNELS HERE ...

shiftTime = [];
if analysisPara.start2zero
    shiftTimeIndx = (startTime ~= 0);
    offset = - mean(startTime(shiftTimeIndx));
    shiftTime = zeros(size(summary));
    for iCML = find(shiftTimeIndx)
        shiftTime(iCML,:) = offset;
    end
end

%% Plot by Calling StandAlone Function

timeCourseAnalysis_StandAlone(summary(goodCML), outputDir, ip.Unmatched, 'shiftTime', shiftTime(goodCML));

%% Save

save([outputDir filesep 'analysisData.mat'], 'directory_CML', 'analysisPara', 'summary');

end

