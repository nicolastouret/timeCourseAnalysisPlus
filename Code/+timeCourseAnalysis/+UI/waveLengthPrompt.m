function [ emissionWL, emissionStr, canceled ] = waveLengthPrompt(cIndx)
%waveLengthPrompt prompt for wavelength selection
%
% OUTPUT
% emissionwL: Number indicating the wavelength selected or empty if canceled
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

% Tae Kim, originally part of UITimeCourseAnalysis
% Adapted into a function by Mark Kittisopikul, Nov 2015

    canceled = false;
    emissionWL = -1;
    emissionStr = [];
    
    emissionStr = inputdlg(['Enter emission wavelength for Channel ' num2str(cIndx)],'Emission wavelength');
    if(isempty(emissionStr))
        canceled = true;
    else
        emissionWL  = str2double(emissionStr);
    end
    
    return;

    stringListWL = { ...
        '525nm : Alexa 488' ,  525
        '530nm : GFP'       ,  530 
        '590nm : Rhod Red X',  590 
        '668nm : Alexa 640' ,  668 
        '669nm : Atto 547N' ,  559 
        'Brightfield'       ,   []
        };
    userChoiceWL = listdlg('PromptString','Select wavelength:', 'SelectionMode','single', 'ListString', stringListWL(:,1));
    
    if(isempty(userChoiceWL))
        % Selection canceled
%         emissionWL = [];
        canceled = true;
    else
        % Selection made;
        emissionWL = stringListWL{userChoiceWL,2};
        if(nargout > 1)
            emissionStr = stringListWL{userChoiceWL,1};
        end
    end
%     if userChoiceWL == 1
%         emissionWL = 525;
%     end
%     if userChoiceWL == 2
%         emissionWL = 530;
%     end
%     if userChoiceWL == 3
%         emissionWL = 590;
%     end
%     if userChoiceWL == 4
%         emissionWL = 668;
%     end
%     if userChoiceWL == 5
%         emissionWL = 669;
%     end
%     if userChoiceWL == 6
%         emissionWL = [];
%         param.imageType_ = '';
%     end
    
end

