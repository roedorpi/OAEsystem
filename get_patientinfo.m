function res = get_patientinfo(pid,folder) 
% pid can be either string or
%   between 1 and 99.
% folder is a string with the 
%   directory to search in.
%    
% res = -1 if there was an error in loading the file.
% res = 0  if the filename was valid but simply could not be found in
%            the specified folder.
% res = HearingHistory struct if the function completed successfully.

    try 
        % pid = '3';                                       % For testing
        % folder = '.';                                    % For testing

        if isempty(pid) == 1; 
            res = -1;
            error('Please provide a patient ID.')
        end
        
        if isnumeric(pid); pid = num2str(pid); end
        if isempty(regexp(pid,'^\d{1,2}$')) == 1
%         if str2num(pid) < 1 || str2num(pid) > 99
            res = -1;
            error('Patient ID not valid.');
        end
        
        if length(pid) == 1
            pid = ['0' pid];
        else%if length(pid) > 2
            %pid = pid(end-1:end);
        end
        extra = '';
        if strcmp(folder(end),'/') ~= 1
            extra = '/';
        end
        folder = [folder extra];
        if exist(folder,'dir') ~= 7
            res = -1;
            error('Specified folder not existing.');
        end
        folder_contents = dir([folder '*.mat']);
        entries_in_fold = length(folder_contents);


        % pid = [char((str2num(pid)<10)*'0') pid];
        
        regex = ['^[FM]' pid '.mat$'];

        % fname = ['F' pid '.mat'];                        % For testing
        found = 0;
        for i = 1:1:entries_in_fold
            fname = folder_contents(i).name;
            if regexp(fname,regex) == 1
                found = 1;
                break;
            end
        end
        if found == 0
            res = -1;
            error('Did not find a patient file with the given ID.')
        end
        
        load([folder fname]) % The HearingHistory struct comes out of this.
        
        res = HearingHistory;% Interesting fields are ID, Name, email, Age,
                             %   Gender, ScreeningPlace, ExposureRating.
    catch err
        res = {res, err.message};
        % display(['get_patientinfo error: ' err.message]);
    end
    
end
