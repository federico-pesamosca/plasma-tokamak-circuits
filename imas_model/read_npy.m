function M = read_npy(filename)
% minimal reader of numpy .npy files (format versions 1.0, 2.0, 3.0)
% dtype '<f8' or '<i8' only, 1 or 2 dimensions; returns a double matrix in
% the stated shape (a 1-D array as a column vector), honouring fortran_order

fid = fopen(filename,'r','ieee-le');
if fid < 0
    error('plasmaless:file','Cannot open %s.', filename);
end
cleanup = onCleanup(@() fclose(fid));

magic = fread(fid,6,'uint8=>double')';
if ~isequal(magic,[147 double('NUMPY')])
    error('plasmaless:file','%s is not a .npy file (bad magic).', filename);
end
ver = fread(fid,2,'uint8=>double');
switch ver(1)
    case 1
        hlen = fread(fid,1,'uint16=>double');
    case {2, 3}
        hlen = fread(fid,1,'uint32=>double');
    otherwise
        error('plasmaless:file','%s: unsupported .npy format version %d.%d.', filename, ver(1), ver(2));
end
header = fread(fid,hlen,'uint8=>char')';

descr = regexp(header,'''descr''\s*:\s*''([^'']*)''','tokens','once');
fortran = regexp(header,'''fortran_order''\s*:\s*(True|False)','tokens','once');
shape = regexp(header,'''shape''\s*:\s*\(([^)]*)\)','tokens','once');
if isempty(descr) || isempty(fortran) || isempty(shape)
    error('plasmaless:file','%s: cannot parse the .npy header ''%s''.', filename, strtrim(header));
end
switch descr{1}
    case '<f8'
        precision = 'float64=>double';
    case '<i8'
        precision = 'int64=>double';
    otherwise
        error('plasmaless:file','%s: unsupported dtype ''%s'' (only <f8 and <i8).', filename, descr{1});
end
dims = sscanf(shape{1},'%d,')';
if numel(dims) < 1 || numel(dims) > 2
    error('plasmaless:file','%s: shape (%s) must have 1 or 2 dimensions.', filename, shape{1});
end

data = fread(fid,prod(dims),precision);
if numel(data) ~= prod(dims)
    error('plasmaless:file','%s: expected %d values, read %d.', filename, prod(dims), numel(data));
end
if numel(dims) == 1
    M = data(:);
elseif strcmp(fortran{1},'True')
    M = reshape(data,dims(1),dims(2));
else
    M = reshape(data,dims(2),dims(1))';
end

end
