function [ job ] = batch( p, cluster , maxWorkers)
%timeCourseAnalysis.run.batch Run timeCourseAnalysis in batch mode
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

    if(nargin < 2 || isempty(cluster))
        cluster = parcluster(p.batchClusterName);
    end
    if(nargin < 3)
        maxWorkers = 24;
    end
    job = batch(cluster,@timeCourseAnalysis, ...
        0, ... % No output
        {p.CML_FullPath, p.outputDir, p}, ...
        'AutoAttachFiles',false, ...
        'CaptureDiary',true, ...
        'Pool',min(cluster.NumWorkers,maxWorkers) - 1 ...
    );

end

