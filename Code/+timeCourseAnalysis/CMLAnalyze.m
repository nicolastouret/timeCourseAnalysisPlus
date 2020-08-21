%deals with individual CML
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
function [CMLSummary, CMLTime, CMLExtra, startTime] = CMLAnalyze(CML,analysisPara)
    alignEvent = CML.analysisPara_.alignEvent;
    %[CMLSummary, CMLTime, CMLExtra] = arrayfun(@(x) MLAnalyze(x, alignEvent), CML.movieLists_, 'UniformOutput', false);
%     nML = numel(CML.movieLists_);
%     CMLSummary = cell(1,nML);
%     CMLTime = cell(1,nML);
%     CMLExtra = cell(1,nML);
%     startTime = zeros(1,nML);
    movieLists = CML.movieLists_;
%     parfor iML = 1:nML
%         [CMLSummary{iML}, CMLTime{iML}, CMLExtra{iML}, startTime(iML)] = timeCourseAnalysis.MLAnalyze(movieLists(iML), alignEvent,analysisPara);
%     end
%     [CMLSummary,CMLTime,CMLExtra,startTime] = ...
%         pararrayfun_progress(...
%         @(x) timeCourseAnalysis.MLAnalyze(x,alignEvent,analysisPara) ...
%         , movieLists ...
%         , 'UniformOutput',false ...
%         , 'DisplayFunc', 'progressTextMultiple' ...
%         );
    [CMLSummary,CMLTime,CMLExtra,startTime] = ...
        arrayfun(...
        @(x) timeCourseAnalysis.MLAnalyze(x,alignEvent,analysisPara) ...
        , movieLists ...
        , 'UniformOutput',false ...
        );
    CMLSummary = vertcat(CMLSummary{:});
    CMLTime = vertcat(CMLTime{:});
    CMLExtra = vertcat(CMLExtra{:});
%     startTime = mean(startTime);
    startTime = mean([startTime{:}]);
end
