% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Creates an unstructured surface-data netCDF file for CLM45.
%
% INPUT:
%       lati_region = Vector containing latitude @ cell-center.
%       long_region = Vector containing longitude @ cell-center.
%       clm_gridded_surfdata_filename = Gridded CLM surface data file
%       out_netcdf_dir = Directory where CLM surface dataset will be saved
%       clm_usrdat_name = User defined name for CLM dataset
%       set_natural_veg_frac_to_one =
%
% Gautam Bisht (gbisht@lbl.gov)
% 01-02-2014
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function fname_out = CreateCLMUgridSurfdatForCLM45(lati_region, long_region, ...
    clm_gridded_surfdata_filename, ...
    out_netcdf_dir, clm_usrdat_name, ...
    set_natural_veg_frac_to_one)

fname_out = sprintf('%s/surfdata_%s_%s.nc',out_netcdf_dir,clm_usrdat_name,datestr(now, 'cyymmdd'));
disp(['  surface_dataset: ' fname_out])

% Check if the file is available
[s,~]=system(['ls ' clm_gridded_surfdata_filename]);

if (s ~= 0)
    error(['File not found: ' clm_gridded_surfdata_filename]);
end

ncid_inp = netcdf.open(clm_gridded_surfdata_filename,'NC_NOWRITE');
ncid_out = netcdf.create(fname_out,'64BIT_OFFSET');

info_inp = ncinfo(clm_gridded_surfdata_filename)

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_inp)

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%                           Define dimensions
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
dimid(1:ndims) = -1;
lonlat_found = 0;
time_found = 0;
in_dim_id = [];
in_dim_name = {};
out_dim_id = [];
out_dim_name = {};

for idim = 1:ndims
    [dimname, dimlen] = netcdf.inqDim(ncid_inp,idim-1);
    disp(['Inp: Dimension name:' dimname])
    in_dim_id = [in_dim_id; idim]
    in_dim_name = [in_dim_name; {dimname}]
    
    switch dimname
        case {'lsmlon','lsmlat'}
            if (lonlat_found == 0) & (dimname == 'lsmlon')
                lonlat_found = 1;
                dimname = 'gridcell';
                dimlen = length(long_region);
                disp(['Out: Dimension name:' dimname])
                dimid(idim) = netcdf.defDim(ncid_out,dimname,dimlen);
                out_dim_name = [out_dim_name; {dimname}];
            end
        case 'gridcell'
            dimlen = length(long_region);
            dimid(idim) = netcdf.defDim(ncid_out,dimname,dimlen);
            out_dim_name = [out_dim_name; {dimname}];
        case 'time'
            time_found = 1;
        otherwise
            disp(['Out: Dimension name:' dimname])
            dimid(idim) = netcdf.defDim(ncid_out,dimname,dimlen);
            out_dim_name = [out_dim_name; {dimname}];
    end
end

if (time_found == 1)
    [time_index] = find(strcmp(in_dim_name,'time'));
    [dimname, dimlen] = netcdf.inqDim(ncid_inp, time_index-1);
    last_dim = ndims + 1
    if(unlimdimid(1) ~= -1)
        nc_time_const = netcdf.getConstant('NC_UNLIMITED')
        dimid(last_dim) = netcdf.defDim(ncid_out,dimname,nc_time_const);
    else
        dimid(last_dim) = netcdf.defDim(ncid_out,dimname,dimlen);
    end
    out_dim_name = [out_dim_name; {dimname}]
end

in_dim_id = num2cell(in_dim_id)
out_dim_id = [1:numel(out_dim_name)]'-1
out_dim_id = num2cell(out_dim_id)
in_dict = containers.Map(in_dim_id, in_dim_name)
out_dict = containers.Map(out_dim_name, out_dim_id)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%                           Define variables
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for ivar = 1:nvars
    [varname,xtype,dim_ids,natts] = netcdf.inqVar(ncid_inp,ivar-1);
    fprintf('\nvarname: %s \ndimids: ', varname);
    fprintf(' %d   ', dim_ids);

    out_dims = map_input_to_output_dimensions(dim_ids, in_dict, out_dict)
    
    varid(ivar) = netcdf.defVar(ncid_out,varname,xtype,out_dims);
    varnames{ivar} = varname;
    %disp([num2str(ivar) ') varname : ' varname ' ' num2str(dimids)])
    
    for iatt = 1:natts
        attname = netcdf.inqAttName(ncid_inp,ivar-1,iatt-1);
        attvalue = netcdf.getAtt(ncid_inp,ivar-1,attname);
        netcdf.putAtt(ncid_out,ivar-1,attname,attvalue);
    end
    
end
varid = netcdf.getConstant('GLOBAL');

[~,user_name]=system('echo $USER');
netcdf.putAtt(ncid_out,varid,'Created_by' ,user_name(1:end-1));
netcdf.putAtt(ncid_out,varid,'Created_on' ,datestr(now,'ddd mmm dd HH:MM:SS yyyy '));
netcdf.endDef(ncid_out);

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Find the nearest neighbor index for (long_region,lati_xy) within global
% dataset
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Get Lat/Lon for global 2D grid.
for ivar = 1:length(varnames)
    if(strcmp(varnames{ivar},'LATIXY'))
        latixy = netcdf.getVar(ncid_inp,ivar-1);
    end
    if(strcmp(varnames{ivar},'LONGXY'))
        longxy = netcdf.getVar(ncid_inp,ivar-1);
    end
end

% read in global pft mask 1=valid 0=invalid
pftmask = ncread(clm_gridded_surfdata_filename,'PFTDATA_MASK');

% mark invalid gridcells as [lon, lat] [-9999, -9999]
latixy(pftmask==0)=-9999;
longxy(pftmask==0)=-9999;

% allocate memoery
ii_idx = zeros(size(long_region));
jj_idx = zeros(size(long_region));

% find the index
for ii=1:size(long_region,1)
    for jj=1:size(long_region,2)
        dist = (longxy - long_region(ii,jj)).^2 + (latixy - lati_region(ii,jj)).^2;
        [nearest_cell_i_idx, nearest_cell_j_idx] = find( dist == min(min(dist)));
        if (length(nearest_cell_i_idx) > 1)
            disp(['  WARNING: Site with (lat,lon) = (' sprintf('%f',lati_region(ii,jj)) ...
                sprintf(',%f',long_region(ii,jj)) ') has more than one cells ' ...
                'that are equidistant.' char(10) ...
                '           Picking the first closest grid cell.']);
            for kk = 1:length(nearest_cell_i_idx)
                disp(sprintf('\t\tPossible grid cells: %f %f', ...
                    latixy(nearest_cell_i_idx(kk),nearest_cell_j_idx(kk)), ...
                    longxy(nearest_cell_i_idx(kk),nearest_cell_j_idx(kk))));
            end
        end
        ii_idx(ii,jj) = nearest_cell_i_idx(1);
        jj_idx(ii,jj) = nearest_cell_j_idx(1);
    end
end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
%                           Copy variables
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for ivar = 1:nvars
    
    %disp(varnames{ivar})
    [varname,vartype,vardimids,varnatts]=netcdf.inqVar(ncid_inp,ivar-1);
    data = netcdf.getVar(ncid_inp,ivar-1);
    varname
    % I believe looking for variable ID == 0 below was meant to find spatial dimensions - this is ambiguous
    % if netcdfs are created across coding languages, or by different groups, the zero dimension will not always be the same
    % I used Python/Xarray, which doesnt fix/garuntee dimension order of the dimension list, to manipulate CLM surface data
    % variable dims are maintained in row- vs column-major dimension order as expected but that is separate from the netcdfs dim list
    % The data variables are also transposed properly by the netcdf API when read into Matlab vs Python as expected
    % but, the dimension list is not fully controllable in python/xarray and can be read out in different orders
    % this is most problematic with the UNLIMITED time dim that is defined first by Xarray but must be defined last in Matlab
    % the work around here has been to define time last and simply map associations between input vs output dim numbers based on dimnames            
    % instead of checking for a dim number of zero, i confirm a spatial dim before converting the data below
    spatial_dims = {'lsmlon', 'lsmlat', 'gridcell'}
    if(isempty(vardimids)==0)
        vdim_names = {};
        for dim_itr = 1:numel(vardimids)
            vdim_id = double(vardimids(dim_itr)+1)
            vdim_names = [vdim_names, in_dict(vdim_id)];
        end
    end
    
    switch varname
        case {'LATIXY'}
            netcdf.putVar(ncid_out,ivar-1,lati_region);
        case {'LONGXY'}
            netcdf.putVar(ncid_out,ivar-1,long_region);
        otherwise
        
            switch length(vardimids)
                case 0
                    netcdf.putVar(ncid_out,ivar-1,data);
                case 1
                    if (lonlat_found)
                        data = 0;
                    else
                        if any(ismember(vdim_names, spatial_dims)) % check if dim is spatial (is this an actual use case? only true for gridcell?)
                            data = data(ii_idx);
                        else
                            data = 0;
                        end
                    end
                    netcdf.putVar(ncid_out,ivar-1,0,length(data),data);
                case 2
                    if any(ismember(vdim_names, spatial_dims)) 
                        
                        if (lonlat_found)
                            data_1d = sgrid_convert_2d_to_1d(vardimids, ii_idx, jj_idx, data, in_dict, out_dict);
                            data_new = PerformFractionCoverCheck(varname, data_1d, set_natural_veg_frac_to_one);
                        else
                            data_2d = ugrid_convert_2d_to_2d(ii_idx, data);
                            data_new = PerformFractionCoverCheck(varname, data_2d, set_natural_veg_frac_to_one);
                        end

                        netcdf.putVar(ncid_out,ivar-1,data_new);
                    else
                        netcdf.putVar(ncid_out,ivar-1,data);
                    end
                case 3
                    if any(ismember(vdim_names, spatial_dims))
                        if (lonlat_found)
                            data_2d  = sgrid_convert_3d_to_2d(vardimids, ii_idx, jj_idx, data, in_dict, out_dict);
                            data_new = PerformFractionCoverCheck(varname, data_2d, set_natural_veg_frac_to_one);
                            netcdf.putVar(ncid_out,ivar-1,data_new);
                        else
                            data_3d  = ugrid_convert_3d_to_3d(ii_idx, data);
                            data_new = PerformFractionCoverCheck(varname, data_3d, set_natural_veg_frac_to_one);
                            netcdf.putVar(ncid_out,ivar-1,zeros(length(size(data_new)),1)',size(data_new),data_3d);
                        end
                        
                    else
                        netcdf.putVar(ncid_out,ivar-1,data);
                    end
                case 4
                    if any(ismember(vdim_names, spatial_dims))
                        if (lonlat_found)
                            data_3d = sgrid_convert_4d_to_3d(vardimids, ii_idx, jj_idx, data, in_dict, out_dict);
                        else
                            disp('error')
                        end
                        
                        netcdf.putVar(ncid_out,ivar-1,zeros(length(size(data_3d)),1)',size(data_3d),data_3d);
                    else
                        netcdf.putVar(ncid_out,ivar-1,data);
                    end
                otherwise
                    disp('error')
            end
    end
end

% close files
netcdf.close(ncid_inp);
netcdf.close(ncid_out);

end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Converts a 2D data (lat,lon) to 1D (gridcell)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function data_1d = sgrid_convert_2d_to_1d(vardimids, ii_idx, jj_idx, data, in_dict, out_dict)

data_2d = zeros(size(ii_idx));

for ii=1:size(ii_idx,1)
    for jj=1:size(jj_idx,2)
        data_2d(ii,jj) = data(ii_idx(ii,jj),jj_idx(ii,jj));
    end
end

% (lon,lat) --> % (gridcell)
%vardimids_new =  [0 vardimids(3:end)-1];
%vardimids = vardimids_new;
vardimids = map_input_to_output_dimensions(vardimids, in_dict, out_dict);
dims = size(data_2d);
if (length(dims)>2)
    dims_new = [dims(1)*dims(2) dims(3:end)];
else
    dims_new = [dims(1)*dims(2) 1];
end
data_1d = reshape(data_2d,dims_new);

end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Converts a 2D data (gridcell,:) to 1D (gridcell,:)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function data_2d = ugrid_convert_2d_to_2d(ii_idx, data)

data_2d = data(ii_idx,:);

end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Converts a 3D data (lat,lon,:) to 2D (gridcell,:)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function data_2d = sgrid_convert_3d_to_2d(vardimids, ii_idx, jj_idx, data, in_dict, out_dict)

nx = size(ii_idx,1);
ny = size(jj_idx,2);
nz = size(data,3);
data_3d = zeros(nx,ny,nz);
for ii = 1:nx
    for jj = 1:ny
        data_3d(ii,jj,:) = data(ii_idx(ii,jj),jj_idx(ii,jj),:);
    end
end

% (lon,lat,:) --> % (gridcell,:)
%vardimids_new =  [0 vardimids(3:end)-1];
%vardimids = vardimids_new;
vardimids = map_input_to_output_dimensions(vardimids, in_dict, out_dict);
dims = size(data_3d);
if (length(dims)>2)
    dims_new = [dims(1)*dims(2) dims(3:end)];
else
    dims_new = [dims(1)*dims(2) 1];
end
data_2d = reshape(data_3d,dims_new);

end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Converts a 3D data (gridcell,:,:) to 3D (gridcell,:,:)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function data_3d = ugrid_convert_3d_to_3d(ii_idx, data)

nx = size(ii_idx,1);
ny = size(data,2);
nz = size(data,3);

data_3d = zeros(nx,ny,nz);
for ii = 1:nx
    data_3d(ii,:,:) = data(ii_idx(ii),:,:);
end

end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Converts a 4D data (lat,lon,:,:) to 3D (gridcell,:,:)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function data_3d = sgrid_convert_4d_to_3d(vardimids, ii_idx, jj_idx, data, in_dict, out_dict)

nx = size(ii_idx,1);
ny = size(ii_idx,2);
nz = size(data,3);
na = size(data,4);

data_4d = zeros(nx,ny,nz,na);
for ii = 1:nx
    for jj = 1:ny
        data_4d(ii,jj,:,:) = data(ii_idx(ii,jj),jj_idx(ii,jj),:,:);
    end
end

% (lon,lat,:) --> % (gridcell,:)
%vardimids_new =  [0 vardimids(3:end)-1];
%vardimids = vardimids_new;
vardimids = map_input_to_output_dimensions(vardimids, in_dict, out_dict);
dims = size(data_4d);
if (length(dims)>2)
    dims_new = [dims(1)*dims(2) dims(3:end)];
else
    dims_new = [dims(1)*dims(2) 1];
end

data_3d = reshape(data_4d,dims_new);

end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Maps between input and output dimension numbers 
% Deals with netcdf surface data with variable dimension order (unlimited time first, etc.)
% in conjunction with changes to the "Define dims" section that simply defines the time dim last
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function out_dims = map_input_to_output_dimensions(dim_ids, in_dict, out_dict)
    
if(isempty(dim_ids)==0)
    vdim_names = {};
    for dim_itr = 1:numel(dim_ids)
        vdim_id = double(dim_ids(dim_itr)+1)
        vdim_names = [vdim_names; in_dict(vdim_id)];
        fprintf('\ndimnames: ');
        fprintf('%s  ', vdim_names{:});
    end
    if any(strcmp(vdim_names,'lsmlon'))
        rm_lonlat = {'lsmlon';'lsmlat'}
        dimupdate = setdiff(vdim_names, rm_lonlat, 'stable')
        if(isempty(dimupdate) == 0)
            diminputs = {};
            diminputs = [diminputs; {'gridcell'}]
            for i = 1:numel(dimupdate)
                diminputs = [diminputs; dimupdate{i}]
            end
        else
            diminputs = {'gridcell'}
        end
    else
        diminputs = vdim_names 
    end
else
    diminputs = [];
end
out_dims = [];
if (isempty(diminputs)==0)
    for dim_itr = 1:numel(diminputs)
        vdim_out = char(diminputs(dim_itr))
        out_dims = [out_dims, out_dict(vdim_out)]
    end
end


end
