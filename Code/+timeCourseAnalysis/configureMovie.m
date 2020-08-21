function [ MD, movieFileName ] = configureMovie( fileName, filePath, param )
%timeCourseAnalysis.configureMovie Configure movie with given parameters
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

name = strsplit(fileName,'.');
name = name{1};

movieFileName = [filePath name filesep name '.mat'];
%create new MD if one doesn't exist already
%use evalc to silence the output
if ~exist(movieFileName, 'file')
    % Movie does not exist
    MD = MovieData([filePath fileName]);
    
    relSize = MD.imSize_./MD.imSize_(1);

    
%KJ 190319: I am commenting out the following lines because they are
%heuristics that were specific to our Spectral microscope images before
%Bioformats could read them. Now they are not necessary for the Spectral
%images because Bioformats can read them properly. More importantly, they
%are causing problems with our Olympus microscope images, because the
%2-channel images from that microscope have this aspect ratio now after
%separating the two channels.

%     if(all(relSize == [ 1 2]))
%         MD = CroppableMovieData.subDivide(MD,relSize,'movieDataFileName_',[name '.mat'],'movieDataPath_',[filePath name]);
%     end
    
    MD.pixelSize_ = param.pixelSize_;
    MD.timeInterval_ = param.timeInterval_;
    MD.numAperture_ = param.numAperture_;
    
    MD.sanityCheck;
    
    for c = 1:length(MD.channels_)
        try
            MD.channels_(c).emissionWavelength_ = param.emissionWavelength_(min(c,end));
        catch err
            warning(['Could not set emissionWavelength_ property to ' num2str(param.emissionWavelength_) ...
                'It is already set to ' MD.channels_(c).emissionWavelength_]);
        end
        MD.channels_(c).exposureTime_ = param.exposureTime_(min(c,end));
        MD.channels_(c).imageType_ = param.imageType_{min(c,end)};
        MD.channels_(c).sanityCheck;
    end
    MD.save;
else
    % Do not override properties above if MovieData already created
    MD = MovieData.load(movieFileName);
    arrayfun(@sanityCheck,MD.channels_,'Unif',false);
end


end

