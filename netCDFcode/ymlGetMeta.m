function metgot = ymlGetMeta(ymlFile,metget)
%
% [metgot] = function ymlGetMeta(ylmFile,metget)
% get gage metadata from .yml file
%
% Inputs:
% ymlFile = .yml file with metadata for this gage
%      e.g. 'D:/Pats_files/FDIF/ocean_templates_fromPike/wind_metadata_d832.yml'
% metget = cell array of metadata values to get, needs to be specific, e.g.
% {'lat','lon','z'}
%
% Output:
%  metgot = values coinciding with metadata specified by metget, e.g.
%  {33.7210, -78.0148, 3.2}
% 
% by Patrick Dickhudt - 25-Feb-2016

fmt = '%s%s';
LF=sprintf('\n');
fid = fopen(ymlFile);
if fid>0
    ind = 1;
    while ~feof(fid)
        nl = fgetl(fid);
        fc = regexp(nl,':');
        C = {nl(1:fc(1)-1),nl(fc(1)+1:end)};
        %C = textscan(nl,fmt,'delimiter',':');
        meta{ind} = C;
        if strfind(meta{ind}{2},' >')
             nl = fgetl(fid);
             newcat = 1;
            while(strcmp(nl(1),' ')) 
                if newcat
                  meta{ind}(2) = {nl};  
                  newcat = 0;
                else
                    meta{ind}(2) = {[meta{ind}{2} LF nl]};
                end
                nl = fgetl(fid);
            end
        end
        ind = ind+1;
    end
    fclose(fid);

    %%
    for i2 = 1:length(metget)
        meti = metget{i2};

        if strcmp(meti,'lat')
            metname = 'geospatial_lat_max';
            fmt = 'd';
        elseif strcmp(meti,'lon')
            metname = 'geospatial_lon_max';
            fmt = 'd';
        elseif strcmp(meti,'Z')
            metname = 'geospatial_vertical_max';
            fmt = 'd';
        elseif strcmp(meti,'serialNumber')
            metname = 'serialNumber';
            fmt = 's';
        elseif strcmp(meti,'pressure_offset')
            metname = 'pressure_offset';
            fmt = 's';
        else
            metname = meti;
            fmt = 'd';
        end

        for i3 = 1:length(meta)
            if strcmp(meta{i3}{1},metname)
                if strcmp(fmt,'d')
                    metgot{i2} = str2num(meta{i3}{2});
                elseif strcmp(fmt,'s')
                    metgot{i2} = meta{i3}{2};
                end
            end
        end
    end
else
    disp('Could not open yml file')
    return
end