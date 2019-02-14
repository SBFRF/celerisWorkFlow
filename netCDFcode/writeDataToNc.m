function writeDataToNc(ncfile, ncid, template, data, dim_map, index, logger)
% This function writes data to netCDF variables
% Inputs:
% 	ncfile   - filename of NetCDF file
%   ncid     - ncid of NetCDF file (ncid = netcdf.create('foo.nc','NC_WRITE');)
%   template - MATLAB structure containing the template defining the NetCDF variables
%	data     - MATLAB structure containing the data. Structure field names should match the variable
% 			   names defined in the template
%   dim_map  - mapping of dimension name and id dim_map.time = dimid
%   index    - index of where to write in NetCDF file (1=create variables)

% % Read template file
% template = ReadYaml(template_file);
%
% Development History
% 04-25-2016    RJF     Adds chunking

if isfield(template, 'x_attributes')
	for n = 1:length(template.x_attributes)
		if isfield(data, template.x_attributes{n})
			% Write some Global Attributes
% 		    varid = netcdf.getConstant('GLOBAL');
		    % netcdf.putAtt(ncfile,varid,template.x_attributes{n},data.(template.x_attributes{n}));
%             netcdf.putAtt(ncid, varid, template.x_attributes{n}, data.(template.x_attributes{n}))
			ncwriteatt(ncfile, '/', template.x_attributes{n}, data.(template.x_attributes{n}))
		end
	end
end

accept_vars = template.x_variables;

for n = 1:length(accept_vars)
    
    % Check if this is a nested structure variable!
    C = strsplit(accept_vars{n}, '.');
    if length(C) == 1
        structure = data;
        var = accept_vars{n};
    else
        structure = data.(C{1});
        var = C{2};
        % For whatever reason the Matlab YAML parser turns the '.' into its HEX ASCII value. 
        % Need to account for that here
        accept_vars{n} = strrep(accept_vars{n}, '.', '0x2E');
    end
    
    if isfield(structure, var)
        var_name = template.(accept_vars{n}).name;
        dimensions = template.(accept_vars{n}).dim;
            
        % If index = 1, create the variable and write the attributes,
        % otherwise just append the data to the open file at index...
        if index == 1
            logger.info('',['Creating ' var '........']);
            createVariable(ncid, template.(accept_vars{n}), var_name, dim_map)
        end

        % Write the data (1D, 2D, or 3D)
        if index == 1 || any(strcmp(dimensions(:),'time'))
            logger.info('',['Writing ' var ' to NetCDF........']);
            ndim = length(dimensions);
            % Check for cell arrays
            if iscell(structure.(var))
                if ndim >= 3 % This must be 2D wave spectra
    %                 for dim = 1:ndim 
    %                     dimid = netcdf.inqDimID(ncid,dimensions{dim});
    %                     [~, dimlen] = netcdf.inqDim(ncid,dimid);
    %                     dims(dim) = dimlen;
    %                 end
    %                 transformed_data = reshape(transformed_data, dims);
                    for i = 1:length(structure.(var))
                        transformed_data(:,:,i) = structure.(var){i}';
                    end
                else % This must be 1D frequency spectra
                    transformed_data = cell2mat(structure.(var)(:))';
                    [r,c] = size(transformed_data);
                    if length(data.(dimensions{1})) ~= c
                        transformed_data = cell2mat(structure.(var)(:)');
                    end
                end
                ncwrite(ncfile, var_name, transformed_data, [ones(1,ndim-1) index]);
                clear transformed_data
          else
               ncwrite(ncfile, var_name, structure.(var), [ones(1,ndim-1) index]);
            end
        end
    else
        logger.info('',[var ' not found in data']);
    end
end


function createVariable(ncid, var, var_name, dim_map)
%--------------------------------------------------------------------------    
% Create the variable
%--------------------------------------------------------------------------
    data_type = getNetCDFDataType(var.data_type);
    dimensions = var.dim;
    dimensions = flip(dimensions);  % Stupid bug in Matlab netCDF library
                                    % that creates dimensions in reverse
                                    % order
    is_time_variable = false;
    if iscell(dimensions) && ismember('time', dimensions) && ~strcmpi(var_name, 'time')
        is_time_variable = true;
    end
    % This list is part of the variable template but wont get written as
    % attributes
    not_var_attr = {'data_type', 'dim', 'fill_value', 'name'};
    
    dimids = [];  % create array of dimensions
    for n = 1:length(dimensions)
        dimids(n) = dim_map.(dimensions{n});
    end
    varid = netcdf.defVar(ncid, var_name, data_type, dimids);
    
    % Add the fill value
    if isfield(var, 'fill_value')
        fill_value = str2double(var.fill_value);
        netcdf.defVarFill(ncid, varid, false, fill_value);
    end
    
    % Add the deflating
    if is_time_variable  % Time variables should be compressed
        netcdf.defVarDeflate(ncid, varid, true, true, 6);
    end

    % Now add the chunking
    chunks = ones(1,length(dimensions));
    for n = 1:length(dimensions)
        dimindex = dim_map.(dimensions{n});
        [name, len] = netcdf.inqDim(ncid, dimindex);
        if len==0
            % Must be an unlimited dimension... Set chunk to large number
            len = 1000;
        end
        if length(dimensions) > 2 && strcmp(name, 'time')
           % Be careful not to set the chunk too high!
           len = 1;
        end
        chunks(n) = len;
    end
    if ~isempty(chunks)
        netcdf.defVarChunking(ncid, varid, 'CHUNKED', chunks);
    end
    
    fnames = fieldnames(var);
    % Write the variable attributes
    for m = 1:length(fnames)
        if ~ismember(fnames{m}, not_var_attr)
            netcdf.putAtt(ncid, varid, fnames{m}, var.(fnames{m}));
%             ncwriteatt(ncfile, var_name, fnames{m}, var.(fnames{m}));
        end
    end

    % Write the units (mandatory)
    netcdf.putAtt(ncid, varid, 'units', var.units);

    % Write the short_name attribute (Mandatory)
    if isfield(var, 'short_name')
        netcdf.putAtt(ncid, varid, 'short_name', var.short_name);
    else
        netcdf.putAtt(ncid, varid, 'short_name', var_name);
    end
    netcdf.endDef(ncid);


function dtype = getNetCDFDataType(data_type)
%--------------------------------------------------------------------------
% This function translates the python netCDF datatypes to Matlab
%--------------------------------------------------------------------------
	switch data_type
		case 'f8'
% 			dtype = 'double';
            dtype = 'NC_DOUBLE';
		case 'f4'
% 			dtype = 'single';
            dtype = 'NC_FLOAT';
		case 'i8'
% 			dtype = 'int64';
            dtype = 'NC_INT64';
		case 'i4'
% 			dtype = 'int32';
            dtype = 'NC_INT';
		case 'i2'
% 			dtype = 'int16';
            dtype = 'NC_SHORT';
		case 'i1'
% 			dtype = 'int8';
            dtype = 'NC_BYTE';
		case 'S1'
% 			dtype = 'char';
            dtype = 'NC_CHAR';
		case 'u8'
% 			dtype = 'uint64';
            dtype = 'NC_UINT64';
		case 'u4'
% 			dtype = 'uint32';
            dtype = 'NC_UINT';
		case 'u2'
% 			dtype = 'uint16';
            dtype = 'NC_USHORT';
		case 'u1'
% 			dtype = 'uint8';
            dtype = 'NC_UBYTE';
		otherwise
% 			dtype = 'double';
            dtype = 'NC_DOUBLE';
	end
