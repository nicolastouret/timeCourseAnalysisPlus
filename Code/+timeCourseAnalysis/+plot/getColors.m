function [ colors ] = getColors( data )
%getColors Get colors for each element in data
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

condColorAll = {...
    [0 0 0],... %black 1
    [1 0 1],... %magenta 7
    [0 0.5 0],... %dark green 3
    [1 0.6 0.3],... %orange 5
    [0 0 1],... %blue 2
    [0 1 0],... %green 6
    [0 1 1],... %cyan 8
    [1 0 0],... %red 4
    [1 1 0],... %yellow 9
    [0.7 0.5 0],... %brown 10
    [0.7 0.7 0.7]... %gray 11
    [0.5 0.5 1],... %purple 12
    [0.3 0.8 1],... %light blue 13
    };
colors = condColorAll(mod(1:numel(data), length(condColorAll) + 1));
colors = reshape(colors,size(data));

end

