function [ncfile, station] = initNcFile(nc_filename, metadata_file)
    % Create the netCDF file and write global attributes

%     dt_now = datestr(now, 29);
    
    % If the file exists, delete it and create a new one
    if exist(nc_filename, 'file') == 2
        delete(nc_filename);
    end
    
    % Create the nc file
    ncfile = netcdf.create(nc_filename, 'NETCDF4');
    
    % Read metadata file
    metadata = ReadYaml(metadata_file);
    
    % Extract the title
    if ~isfield(metadata, 'title')
        error('title is missing from metadata template!')
    else
        station = metadata.title;
    end
 
    % Get all the metadata fields
    fields = fieldnames(metadata);
    
    % And write the metadata as global attributes
    varid = netcdf.getConstant('GLOBAL');
    for i = 1:numel(fields)
        netcdf.putAtt(ncfile, varid, fields{i}, metadata.(fields{i}));
%         disp(metadata.(fields{i}))
    end
    