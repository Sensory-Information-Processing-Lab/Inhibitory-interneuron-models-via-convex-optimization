function write_log_file(log_file, time_stamp, pars, varargin)
% Write log file
% Inputs:
%   pars: cell array of parameter structs;
%   varargin:
%     'sort_pars': whether to sort the parameter field names
%     alphabetically.
%     'write_parnames': whether to output extra line of  par names to log

% Mengchen Zhu

%% Parse inputs
% Default parameter values
% 0 initial conidtions
sort_pars = false;
write_parnames = false;
options = struct('sort_pars', sort_pars, 'write_parnames', write_parnames);
options = parse_inputs(options, varargin);
% $$$ 
% $$$ % read the acceptable names
% $$$ optionNames = fieldnames(options);
% $$$ 
% $$$ % count arguments
% $$$ nArgs = length(varargin);
% $$$ if round(nArgs/2)~=nArgs/2
% $$$    error('needs propertyName/propertyValue pairs')
% $$$ end
% $$$ 
% $$$ for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
% $$$    inpName = lower(pair{1}); %# make case insensitive
% $$$ 
% $$$    if any(strmatch(inpName,optionNames))
% $$$       options.(inpName) = pair{2};
% $$$    else
% $$$       error('%s is not a recognized parameter name',inpName)
% $$$    end
% $$$ end

% The parameter names
par_name = {'time_stamp'};
% Parameter values
par_values = {time_stamp};
pars_sorted = [];
for idx = 1:length(pars)
    if options.sort_pars
        % Sort the fields in the structure
        pars_sorted = orderfields(pars{idx});
    else
        pars_sorted = pars{idx};
    end
    % The parameter names
    par_name = [par_name, fieldnames(pars_sorted)'];
    par_values = [par_values, struct2cell(pars_sorted)'];
end

% Write log file
if(~exist(log_file, 'file'))
    % If not, generate the file and the first row (parameter
    % names): strip out the field names as column headings
    log_cell = par_name;
    cell2csv(log_file, log_cell);
else
    % If exists, open it and read to a cell.
    log_cell = csv2cell(log_file, 'fromfile'); 
    % If the par_name does not agree with the old one, add a new
    % line of parameter names and new dimensions.
    if (size(log_cell, 2)<size(par_name,2)) | options.write_parnames
        log_cell = [log_cell, cell(size(log_cell,1), (size(par_name,2) ...
                                                      - size(log_cell,2))); ...
                    par_name];
    end
end
% Write the time stamp and the parameters to the last row. 
log_cell = [log_cell; par_values];
% Write the new log_file
cell2csv(log_file, log_cell);
