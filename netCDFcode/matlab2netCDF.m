function matlab2netCDF(data, metadata, template, index, output_file)
% Converts a wavedat MATLAB structure to a NetCDF file

% Inputs:
% 	data - MATLAB structure containing all data needed
%   metadata - Path to YAML file containing the global metadata
% 	template - Path to YAML file containing the template used to define NetCDF schema
%   index - Specifies where to append the data (1 = create NetCDF file)
%	output_file - Absolute path to save netCDF file
%
% 2015-08-27 RJF Updated to work with nested structures
%%
% forcing to a character array if came in as string
metadata = char(metadata);
template = char(template); 

% Setup Logger
if ~(exist(fullfile(pwd,'log'), 'dir') == 7)
    mkdir(fullfile(pwd,'log'))
end
log_file = fullfile(pwd, 'log', 'netcdf_conversion.log');
logger = log4m.getLogger(log_file);
logger.info('','Converting to NetCDF');

% Set up YAMLMatlab toolbox
if ~exist('ReadYaml.m', 'file')
    if exist('YAMLMatlab_0.4.3', 'dir')
        addpath(genpath('YAMLMatlab_0.4.3'))
    else
        error('ReadYaml.m not found!  Please add YAMLMatlab to your search path');
    end
end

% Check inputs
[pathstr, ~, ext] = fileparts(output_file);
if ~strcmp(ext, '.nc')
    error('output_file is a NetCDF file and needs extension .nc')
end
if ~isempty(pathstr) && ~exist(pathstr, 'dir')
    mkdir(pathstr)
end

try
	% Read template file
	nc_template = ReadYaml(template);

	% Check the timezone (need it in UTC)
    if isfield(data, 'timezone')
		tz = data.timezone;
	else
		tz = 0;
    end

	% Get the filename
    % nc_filename = create_nc_filename(filename);
    nc_filename = output_file;

%     % Get the keywords
%     keywords = nc_template.x_keywords;

    % Create the netCDF file and add global attributes
    if index == 1
        logger.info('',['Creating file ' nc_filename]);
        [ncid, station] = initNcFile(nc_filename, metadata);
        % Add a couple extra global variables
        today = datestr(now,29);
        varid = netcdf.getConstant('GLOBAL');
        netcdf.putAtt(ncid, varid, 'date_created', today);
        netcdf.putAtt(ncid, varid, 'date_issued', today);
        netcdf.putAtt(ncid, varid, 'time_coverage_start', datestr(data.time(1),'yyyy-mm-ddTHH:MM:SS'));
        netcdf.putAtt(ncid, varid, 'time_coverage_end', datestr(data.time(end),'yyyy-mm-ddTHH:MM:SS'));
        
        % Add station_name to wavedat
        data.station_name = station;

        % Create the dimensions
        dims = nc_template.x_dimensions;
        
        % Check if this is realtime file
        try
            processing = netcdf.getAtt(ncid, varid, 'processing');
            if strncmpi(processing, 'historic', 6)
                realtime = false;
            else
                realtime = true;
            end
        catch
            realtime = true;
        end
        for n = 1:length(dims)

            logger.info('',['Creating Dimension ' dims{n}]);
            % A couple special cases for dimensions time and station_name_length
            
            if strcmpi(dims(n), 'time') && realtime
                dim_name = 'time';
                dimid = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
            elseif strcmpi(dims(n), 'station_name_length')
                dim_name = 'station_name_length';
                dimid = netcdf.defDim(ncid,'station_name_length',length(data.station_name));
            else
%                 all_var_names = fieldnames(data);
%                 if ~any(ismember(all_var_names, dims{n}))
%                     error('Dimension name is not a variable in data!')
%                 end

                % Check if this is a nested structure variable!
                C = strsplit(dims{n}, '.');
                if length(C) == 1
                    structure = data;
                    var = dims{n};
                else
                    structure = data.(C{1});
                    var = C{2};
                    % For whatever reason the Matlab YAML parser turns the '.' into its HEX ASCII value. 
                    % Need to account for that here
                    dims{n} = strrep(dims{n}, '.', '0x2E');
                end
                if isfield(structure, var)
                    dim_name =nc_template.(dims{n}).name;
                    dim_size = length(structure.(var));
                    dimid = netcdf.defDim(ncid, dim_name, dim_size);
                else
                    logger.info('',['Skipping Dimension ' dims{n} ' - Not found in data']);
                end
            end
            dim_mapping.(dim_name) = dimid;
        end
        netcdf.sync(ncid);
    else
        dim_mapping = struct;  % Make an empty structure
        ncid = netcdf.open(nc_filename,'WRITE');
    end
    
    % Write the data to the netCDF file
    logger.info('',['Writing to file ' nc_filename]);
    
    % Convert time stamp from matlab to UNIX
    data.time = convertToUnixEpoch(data.time, tz);
    
    if isfield(data,'time_QA')
        data.time_QA = convertToUnixEpoch(data.time_QA, tz);
    end
    
    writeDataToNc(nc_filename, ncid, nc_template, data, dim_mapping, index, logger);
    
    % Close the netCDF file
    netcdf.close(ncid);
    
    logger.info('','Write operation successful');

catch ME
    logger.info('','Write operation failed!!!');
    % Log the error
    error_msg = getReport(ME,'extended','hyperlinks','off');
    logger.error('waves2netCDF', error_msg)

    if exist('ncid', 'var')
        % Close the netCDF file
        netcdf.close(ncid);
    end

end