function options = parse_inputs(options, varargin)

% read the acceptable names
optionNames = fieldnames(options);

% count arguments
nArgs = length(varargin{1});
if round(nArgs/2)~=nArgs/2
   error('needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin{1},2,[]) %# pair is {propName;propValue}
   inpName = pair{1}; % case sensitive

   if any(strmatch(inpName,optionNames))
      options.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end
