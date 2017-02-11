function varargout = VestBMS_defaults(command,varargin)
%PROJECT_DEFAULTS Get default variables for a given project.
%
%   OPTLIST = PROJECT_DEFAULTS('options') returns a project-dependent 
%   options list. OPTLIST is a struct array with the following fields:
%      'name'       name of the option
%      'type'       'flag','matrix','string' (flags can take values 0 or 1)
%      'default'    default value
%
%   [MODELSTRING,DATAIDSTRING] = PROJECT_DEFAULTS('strings',MODEL,DATAID) 
%   returns the project-dependent model string and data id string for
%   a given MODEL and DATAID.
%
%   See also MODELWORK_DEFAULTS, PARSEOPTS.

switch lower(command)
    case 'options'
        NTRIMTRIALS = 0; % Trim this number of first trials (per condition)
        NTYPES = 8; % Number of different prior distributions
        NBLOCKSPERSESSION = 8;

        optlist(1) = struct('name', 'nc', 'type', 'matrix', 'default', 3);
        optlist(end+1) = struct('name', 'ntypes', 'type', 'matrix', 'default', NTYPES);
        optlist(end+1) = struct('name', 'runtimemax', 'type', 'matrix', 'default', Inf);        
        optlist(end+1) = struct('name', 'blocksflag', 'type', 'matrix', 'default', true(1, NBLOCKSPERSESSION));
        optlist(end+1) = struct('name', 'blockflags', 'type', 'matrix', 'default', []);
        optlist(end+1) = struct('name', 'trimfirstntrials', 'type', 'matrix', 'default', NTRIMTRIALS);
        optlist(end+1) = struct('name', 'trialtypeflags', 'type', 'matrix', 'default', true(1, NTYPES));
        optlist(end+1) = struct('name', 'priorblocks', 'type', 'flag', 'default', 1);
        optlist(end+1) = struct('name', 'nominimize', 'type', 'flag', 'default', 0);
        optlist(end+1) = struct('name', 'discardcategorical', 'type', 'flag', 'default', 1);
        optlist(end+1) = struct('name', 'discardnoncategorical', 'type', 'flag', 'default', 0);
        optlist(end+1) = struct('name', 'discardunity', 'type', 'flag', 'default', 0);
        % optlist(end+1) = struct('name', 'keeponlyhalf', 'type', 'matrix', 'default', 0);        
        optlist(end+1) = struct('name', 'sessions', 'type', 'matrix', 'default', []);
        optlist(end+1) = struct('name', 'discardlargedisparities', 'type', 'flag', 'default', 0);
        optlist(end+1) = struct('name', 'logpriorflag', 'type', 'flag', 'default', 1);
        % optlist(end+1) = struct('name', 'debugtheta', 'type', 'matrix', 'default', []);
        optlist(end+1) = struct('name', 'unifulltheta', 'type', 'string', 'default', []);
        optlist(end+1) = struct('name', 'unifullthetanumber', 'type', 'matrix', 'default', 1);
        optlist(end+1) = struct('name', 'approxflag', 'type', 'flag', 'default', 0);
        optlist(end+1) = struct('name', 'fakedatastartx', 'type', 'flag', 'default', 0);
        optlist(end+1) = struct('name', 'bincenters', 'type', 'matrix', 'default', []);
        optlist(end+1) = struct('name', 'bincenters_bim', 'type', 'matrix', 'default', []);
        optlist(end+1) = struct('name', 'respbincenters', 'type', 'matrix', 'default', []);
        optlist(end+1) = struct('name', 'binweights', 'type', 'matrix', 'default', []);
        optlist(end+1) = struct('name', 'defaultprior', 'type', 'matrix', 'default', []);
        optlist(end+1) = struct('name', 'loadinitfromconst', 'type', 'matrix', 'default', 0);
        optlist(end+1) = struct('name', 'loadinitfromdisc', 'type', 'matrix', 'default', 0);
        
        varargout{1} = optlist;
        
    case 'strings'        
        model = varargin{1};
        dataid = varargin{2};
        
        modelstring = packuint(model);
        suffix = packuint(dataid(2:end));
        if isempty(suffix)
            dataidstring = ['S' num2str(dataid(1))];                
        else
            dataidstring = ['S' num2str(dataid(1)) '-' packuint(dataid(2:end))];
        end
        
        varargout{1} = modelstring;
        varargout{2} = dataidstring;
        
    case 'plots'
        % plots.NoiseColors = [1 0 0; 0 1 0; 0 0 1];
        
        % plots.NoiseColors = [166,206,227; 31,120,180; 178,223,138]*0.8/255;
        % plots.NoiseColors = [27,158,119;  117,112,179; 217,95,2 ]/255;
        plots.NoiseColors = [66 39 0; 178 130 24; 255 210 56]/255;

        
        varargout{1} = plots;
end