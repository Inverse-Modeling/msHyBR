function val = sshybr_get(options,name,default,flag)
%
%   VAL = sshybr_get(options,name,default,flag)
%
%   sshybr_get gets options parameters.
%   VAL = sshybr_get(OPTIONS,'NAME') extracts the value of the named parameter
%   from sshybr options structure OPTIONS, returning an empty matrix if
%   the parameter value is not specified in OPTIONS.  It is sufficient to
%   type only the leading characters that uniquely identify the
%   parameter.  Case is ignored for parameter names.  [] is a valid OPTIONS
%   argument.
%
%   VAL = sshybr_get(OPTIONS,'NAME',DEFAULT) extracts the named parameter as
%   above, but returns DEFAULT if the named parameter is not specified
%   in OPTIONS.  For example
%
%     param = sshybr_get(opts,'RegPar','WGCV');
%
%   returns param = 'WGCV' if the RegPar property is not specified in opts.
%
%
%   J.Chung 12/2021
%   Modified my M. Sabate Landman 04/2024

% undocumented usage for fast access with no error checking
if (nargin == 4) && isequal(flag,'fast')
    val = sshybr_getfast(options,name,default);
    return
end

if nargin < 2
    error('Not enough input arguments.');
end
if nargin < 3
    default = [];
end

if ~isempty(options) && ~isa(options,'struct')
    error('First argument must be an options structure created with sshybr_set.');
end

if isempty(options)
    val = default;
    return;
end
allfields = {'InSolv'; 'RegPar';'nLevel';'Omega';'Iter';'Reorth'; ...
    'x_true';'BegReg'; 'FlatTol'; 'MinTol'; 'ResTol'; 'thr'; 'mask'};

Names = allfields;

name = deblank(name(:)'); % force this to be a row vector
j = find(strncmpi(name,Names,length(name)));
if isempty(j)               % if no matches
    error(['Unrecognized property name ''%s''.  ' ...
        'See sshybr_set for possibilities.'], name);
elseif length(j) > 1            % if more than one match
    % Check for any exact matches (in case any names are subsets of others)
    k = find(strcmpi(name,Names));
    if length(k) == 1
        j = k;
    else
        msg = sprintf('Ambiguous property name ''%s'' ', name);
        msg = [msg '(' Names{j(1),:}];
        for k = j(2:length(j))'
            msg = [msg ', ' Names{k,:}];
        end
        msg = sprintf('%s).', msg);
        error(msg);
    end
end

if any(strcmp(Names,Names{j,:}))
    val = options.(Names{j,:});
    if isempty(val)
        val = default;
    end
else
    val = default;
end

%------------------------------------------------------------------
function value = sshybr_getfast(options,name,defaultopt)
%HYBRGETFAST- Get HyBR_lsmr OPTIONS parameter with no error checking.
%   VAL = sshybr_getFAST(OPTIONS,FIELDNAME,DEFAULTOPTIONS) will get the
%   value of the FIELDNAME from OPTIONS with no error checking or
%   fieldname completion. If the value is [], it gets the value of the
%   FIELDNAME from DEFAULTOPTIONS, another OPTIONS structure which is
%   probably a subset of the options in OPTIONS.
%
if isempty(options)
     value = defaultopt.(name);
     return;
end
% We need to know if name is a valid field of options, but it is faster to use 
% a try-catch than to test if the field exists and if the field name is
% correct.
try
    value = options.(name);
catch
    value = [];
    lasterr('');  % clean up last error
end

if isempty(value)
    value = defaultopt.(name);
end


