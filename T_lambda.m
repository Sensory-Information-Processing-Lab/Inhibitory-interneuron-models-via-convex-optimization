function a = T_lambda(u, pars, varargin)
% Thresholding function for lca; see Jul 3, 2011 notebook.
% u: state variable
% pars:
%   threshold: a POSITIVE threshold
%   thr_func: type of thresholding function
%   basis_size: normal or double basis
% varargin:
%   'u_block': l2 norms of the block when calculating block
%     l1. Same length as u.

%% Parse inputs
% Default parameter values
u_block = NaN;
options = struct('u_block', u_block);

% read the acceptable names
optionNames = fieldnames(options);

% count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
   inpName = lower(pair{1}); %# make case insensitive

   if any(strmatch(inpName,optionNames))
      options.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end

%% Compute
switch pars.thr_func
  case 'soft'
    switch pars.basis_size
      case 'double'
        % Use one-sided thresholding
        a = u-pars.threshold;
        a(u < pars.threshold) = 0;
      case 'normal'
        % Use two-sided thresholding
        a = u-pars.threshold;
        a(-pars.threshold < u & u < pars.threshold) = 0;
        a(u <= -pars.threshold) = a(u <= -pars.threshold) + 2*pars.threshold;
    end
  case 'hard'
    switch pars.basis_size
      case 'double'
        % Use one-sided thresholding
        a = u;
        a(u < pars.threshold) = 0;
      case 'normal'
        % Use two-sided thresholding
        a = u;
        a(-pars.threshold < u & u < pars.threshold) = 0;        
    end
  case 'block'
    % Block (a.k.a group) l1
    switch pars.basis_size
      case 'normal'
        if(isnan(options.u_block))
            error('Need to provide block l2 norm.');
        else
            a = u - pars.threshold * u ./ options.u_block;
            a(options.u_block < pars.threshold) = 0;
        end
    end
end
