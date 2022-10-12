function [status, fileHead, folderName, dataFileName] = save_result(case_id, string_remark, string_variable_to_save)
% Save the source file (simulation script) and the data file (simulation
% results).
%
% save_result
% save_result(string_variable_to_save)
% save_result(case_id)
% save_result(case_id, string_remark)
% save_result(case_id, [], string_variable_to_save)
% save_result(case_id, string_remark, string_variable_to_save)
%
% The data file is named as {date}-{time}-{remark_string}.mat, and stored
% under ./Results/ by default. The source file is copied to ./Results/ by
% default and renamed as '{original_name}-{date}-{time}-{remark_string}.m'.
%
% USAGE EXAMPLES
% ==============
% save_result(string_variable_to_save) Specifies the variables to save. If 
%           not specify, save all variables.
%
% save_result(1) Saves the simulation results under folder ./Results/Case01/
%
% save_result(2, 'test') Saves the simulation results under folder 
%           ./Results/Case02/ with the data file name being 
%           '{date}-{time}-test.mat'.
%
% save_result(-1, 'test') Saves the simulation results under folder 
%           ./Results/ with the data file name being 
%           '{date}-{time}-test.mat'.

status = 0;

switch nargin 
    case 0
        folderName = 'Results/';
        string_remark = [];
        string_variable_to_save = [];
    case 1
        if isnumeric(case_id)
            folderName = sprintf('Results/Case%02d', case_id);
            string_remark = [];
            string_variable_to_save = [];
        else
            folderName = 'Results/';
            string_variable_to_save = case_id;
            case_id = -1;
            string_remark = [];
        end
            
    case 2
        % warning('The case with only 2 arguments is not recommended (hidden ambiguous)! The string will be treated as remark.');
        folderName = sprintf('Results/Case%02d', case_id);
        string_variable_to_save = [];
        
    case 3
        if case_id >= 0
            folderName = sprintf('Results/Case%02d', case_id);
        else
            folderName = 'Results/';
        end
end
if ~isempty(string_remark)
    string_remark = ['-' string_remark];
    string_remark = regexprep( string_remark, ' ', '_' );
end

if exist(folderName, 'dir') ~= 7
    mkdir(folderName);
end

% Create file name
currentTime = datestr(now, 31);
currentTime = regexprep( currentTime, ' ', '_' );
currentTime = regexprep( currentTime, ':', '-' );

dataFileName = sprintf('%s/%s%s.mat', folderName, currentTime, string_remark);
fileHead = sprintf('%s/%s', folderName, currentTime);

cnt = 0;
MAXCNT = 100;
while exist(dataFileName, 'file') == 2 && cnt < MAXCNT
    cnt = cnt + 1;
    dataFileName = sprintf('%s/%s%s(%d).mat', folderName, currentTime, ...
        string_remark, cnt);
end
if cnt >= MAXCNT
    warning('Failed to create a valid data file!');
    status = -1;
end

% Save data file
save_expression = sprintf('save %s %s', dataFileName, string_variable_to_save);
evalin('caller', save_expression);

% Save source file
[st_caller, ~] = dbstack('-completenames');
sourceFile = st_caller(2).file;
sourceFileName = st_caller(2).name;
sourceFileNewPath = sprintf('%s/%s%s_%s.m', folderName, ...
    currentTime, string_remark, sourceFileName);
copyfile(sourceFile, sourceFileNewPath);
