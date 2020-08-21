function [ job ] = createJob( p, cluster , maxWorkers )
%timeCourseAnalysis.run.createJob Create a parallel.Job for
%timeCourseAnalysis
%
% See also timeCourseAnalysis.run.batch
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
    numWorkers = min(cluster.NumWorkers,maxWorkers);
    job = createCommunicatingJob( ...
        cluster, ...
        'Type','pool', ...
        'NumWorkersRange',[min(6,numWorkers) numWorkers], ...
        'AutoAttachFiles',false, ...
        'Name','Time Course Analysis', ...
        'Tag','UITimeCourseAnalysis' ...
        );
    createTask(job,@timeCourseAnalysis, ...
        0, ... % No output
        {p.CML_FullPath, p.outputDir, p}, ...
        'CaptureDiary', true ...
        );
    
    % Need to submit(job)
    % Get diary by job.Tasks(1).Diary

end

