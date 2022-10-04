function pars_str = pars2str(pars)
% Convert a parameter variable or a parameter structure into a string
% of the following format: field1_val1_field2_val2...  This is useful
% for appending parameter values used in a simulation to file names,
% when the parameters are stored in a struct.

% Accepted input pars field includes: string, numerals, and numeric
% vector.

% Mengchen Zhu


if isstruct(pars)
    % Extract the field names and values
    names = fieldnames(pars);
    val =  struct2cell(pars);
else
    names = {inputname(1)};
    val = {pars};
end
% Interlace field names with their values
pars_str = '';
for n = 1:length(val)
    % convert to string
    if isnumeric(val{n})
        if (length(val{n}) == 1)
            pars_str = [pars_str, names{n}, '_', num2str(val{n})];
        else 
            % numeric array, store the first, the increment, and the last value.
            pars_str = [pars_str, names{n}, '_', num2str(val{n}(1)), ...
                        '_', num2str(val{n}(2) - val{n}(1)), '_', ...
                        num2str(val{n}(end))];
        end        
    elseif islogical(val{n})
        switch val{n}
          case true
            val{n} = 'true';
          case false
            val{n} = 'false';
        end
        pars_str = [pars_str, names{n}, '_', val{n}];
    else
        pars_str = [pars_str, names{n}, '_', val{n}];
    end    
    % Replace "." with "p" to avoid file accessing issue
    pars_str = regexprep(pars_str, '\.','p');
    % No trailing '_' if the last parameter
    if n < length(val)
        pars_str = [pars_str, '_'];
    end
end

