function [model, x0] = configure_linear_model(model_dir, dt, pf_active_md)
% build the discrete plasmaless model from an external continuous linear model
%   dI/dt = A I + B U,   Y = C I
% (e.g. a NICE linearisation) described by linear_model.json in model_dir.
% I (states), U (model inputs) and Y (observables) are deviations from the
% linearisation point; the absolute values are model.op + deviation.
% Descriptor fields (README.md, "External linear model"): description, shot,
% time, files {A, B, C: basenames, .npy or whitespace .txt}, states (n names),
% inputs (m input names, column order of B; e.g. generator voltages),
% input_map (input -> {coil short name: coefficient}: the input as a linear
% combination of the pf_active coil voltages; an input absent or with an
% empty object is not drivable through pf_active: zero row, warning
% plasmaless:inputUndriven; without input_map, inputs must be the coil set
% and the map is the permutation to machine-description order), coil_states
% (coil -> state), passive_states (list), outputs (p names in {rmag, zmag,
% rX, zX, Ip, gap<k>, Dist2ndSep}), optional operating_point {coil_current,
% voltage, outputs} (coil space, in this wrapper's (IMAS) convention), optional
% convention {state_sign, input_sign, output_sign {output: sign}} (each +1 or
% -1, default +1), notes.
% Sign convention: the matrix files are written in the provider's frame,
%   x_file = state_sign * x,  u_file = input_sign * u,  y_file(k) = output_sign(k) * y(k)
% with x, u, y in the wrapper's (IMAS) convention. Substituting into
% dx_file/dt = A x_file + B u_file, y_file = C x_file gives the wrapper-frame
% matrices, applied before the input map and the discretisation:
%   A unchanged,  B_conv = (input_sign/state_sign) * B,
%   C_conv(k,:) = (state_sign/output_sign(k)) * C(k,:).
% (Received model: state_sign = input_sign = -1, output_sign(Ip) = -1, geometric
% outputs +1: A and B unchanged, the Ip row unchanged, the geometric rows negated.)
% The raw file matrices are kept in model.A_file, B_file, C_file.
% pf_active_md: machine-description pf_active; its coil names/order define
% model.coil_names; the model input matrix in coil space is B_coil = B*G with
% G (m x n_coils) built from input_map, and the discrete system is c2d of
% (A, B_coil, C): the step function takes the 14 coil voltages.
% Output model: same fields as configure_plasmaless_model (kind 'linear',
% sysd, A, B, C, D, dt, n_states, state_names, coil_names, n_coils, x_to_coil,
% loop_names, n_loops, x_to_loop, y_names, C_y, op) plus input_names,
% input_map (G), source (model_dir, files, shot, time, description), A_c,
% B_c (continuous, coil space), B_c_input (continuous, input order), C_c (all
% in the wrapper's convention), convention (state_sign, input_sign, output_sign
% p x 1), A_file, B_file, C_file (as read, provider's frame).
% x0: zeros (the run starts at the linearisation point).

if nargin ~= 3
    error('plasmaless:nargin','Use configure_linear_model(model_dir, dt, pf_active_md).');
end
desc_file = fullfile(model_dir,'linear_model.json');
if ~isfile(desc_file)
    error('plasmaless:descriptor','Descriptor %s not found.', desc_file);
end
d = jsondecode(fileread(desc_file));
required = {'description','shot','time','files','states','inputs','coil_states', ...
    'passive_states','outputs'};
for ii = 1:numel(required)
    if ~isfield(d,required{ii})
        error('plasmaless:descriptor','%s: field ''%s'' missing.', desc_file, required{ii});
    end
end
for f = {'A','B','C'}
    if ~isfield(d.files,f{1})
        error('plasmaless:descriptor','%s: files.%s missing.', desc_file, f{1});
    end
end
states = name_list(d.states,'states',desc_file);
inputs = name_list(d.inputs,'inputs',desc_file);
passive = name_list(d.passive_states,'passive_states',desc_file);
outputs = name_list(d.outputs,'outputs',desc_file);
n = numel(states); m = numel(inputs); p = numel(outputs);
if ~isstruct(d.coil_states)
    error('plasmaless:descriptor','%s: coil_states must be an object coil -> state.', desc_file);
end

% outputs vocabulary (before the convention block, which refers to the outputs)
for ii = 1:p
    if isempty(regexp(outputs{ii},'^(rmag|zmag|rX|zX|Ip|Dist2ndSep|gap[1-9][0-9]*)$','once'))
        error('plasmaless:observable', ['%s: output ''%s'' is not in the vocabulary ', ...
            '{rmag, zmag, rX, zX, Ip, gap<k>, Dist2ndSep}.'], desc_file, outputs{ii});
    end
end

% matrices
A = load_matrix(fullfile(model_dir,d.files.A));
B = load_matrix(fullfile(model_dir,d.files.B));
C = load_matrix(fullfile(model_dir,d.files.C));
if ~isequal(size(A),[n n]) || ~isequal(size(B),[n m]) || ~isequal(size(C),[p n])
    error('plasmaless:size', ['%s: A %s, B %s, C %s do not match %d states, %d inputs, ', ...
        '%d outputs.'], desc_file, mat2str(size(A)), mat2str(size(B)), mat2str(size(C)), n, m, p);
end
if ~all(isfinite(A(:))) || ~all(isfinite(B(:))) || ~all(isfinite(C(:)))
    error('plasmaless:size','%s: the matrices contain non-finite values.', desc_file);
end

% sign convention of the files relative to the wrapper (IMAS) frame (header)
conv = sign_convention(d, outputs, desc_file);
A_file = A; B_file = B; C_file = C;
B = (conv.input_sign/conv.state_sign)*B;
C = (conv.state_sign./conv.output_sign).*C;   % row k scaled by state_sign/output_sign(k)

% coils: machine-description names and order
n_coils = numel(pf_active_md.coil);
coil_names = cell(1,n_coils);
for ii = 1:n_coils
    coil_names{ii} = plasmaless_coil_key(pf_active_md.coil{ii});
end
if numel(unique(coil_names)) ~= n_coils
    error('plasmaless:coilName','Duplicate pf_active coil names.');
end
coil_keys = fieldnames(d.coil_states)';
if ~isequal(sort(coil_keys),sort(coil_names))
    error('plasmaless:descriptor','%s: coil_states keys {%s} are not the machine-description coils {%s}.', ...
        desc_file, strjoin(coil_keys,','), strjoin(coil_names,','));
end

% input map G (m x n_coils): input k = sum_j G(k,j) * coil voltage j
G = input_matrix(d, inputs, coil_names, desc_file);
B_coil = B*G;

% states: every state is a coil state or a passive state, not both
x_to_coil = zeros(n_coils,n);
is_coil_state = false(1,n);
for ii = 1:n_coils
    s = d.coil_states.(coil_names{ii});
    if ~ischar(s) && ~isstring(s)
        error('plasmaless:descriptor','%s: coil_states.%s must be a state name.', desc_file, coil_names{ii});
    end
    k = find(strcmp(states,char(s)),1);
    if isempty(k)
        error('plasmaless:descriptor','%s: coil_states.%s = ''%s'' is not a state.', ...
            desc_file, coil_names{ii}, char(s));
    end
    x_to_coil(ii,k) = 1;
    is_coil_state(k) = true;
end
is_passive = ismember(states,passive);
missing = setdiff(passive,states);
if ~isempty(missing)
    error('plasmaless:descriptor','%s: passive_states {%s} are not states.', desc_file, strjoin(missing,','));
end
both = states(is_coil_state & is_passive);
if ~isempty(both)
    error('plasmaless:descriptor','%s: states {%s} are both coil and passive states.', desc_file, strjoin(both,','));
end
none = states(~is_coil_state & ~is_passive);
if ~isempty(none)
    error('plasmaless:descriptor','%s: states {%s} are neither coil nor passive states.', desc_file, strjoin(none,','));
end
loop_names = states(is_passive);
n_loops = numel(loop_names);
x_to_loop = zeros(n_loops,n);
x_to_loop(:,is_passive) = eye(n_loops);

% discrete model (zero-order hold, as configure_plasmaless_model)
ev = eig(A);
if any(real(ev) > 1e-9)
    warning('plasmaless:unstable','%s: %d eigenvalue(s) of A with real part > 1e-9 (max %.3g).', ...
        desc_file, nnz(real(ev) > 1e-9), max(real(ev)));
end
sys = ss(A,B_coil,C,zeros(p,n_coils));
sysd = c2d(sys,dt);
if any(sysd.D(:) ~= 0)
    error('plasmaless:size','%s: the discrete D matrix is not zero.', desc_file);
end

model.kind = 'linear';
model.sysd = sysd;
model.A = sysd.A;
model.B = sysd.B;
model.C = sysd.C;
model.D = sysd.D;
model.dt = dt;
model.n_states = n;
model.state_names = states;
model.coil_names = coil_names;
model.n_coils = n_coils;
model.x_to_coil = x_to_coil;
model.loop_names = loop_names;
model.n_loops = n_loops;
model.x_to_loop = x_to_loop;
model.y_names = outputs;
model.C_y = sysd.C;
model.input_names = inputs;
model.input_map = G;
model.A_c = A;
model.B_c = B_coil;
model.B_c_input = B;
model.C_c = C;
model.convention = conv;
model.A_file = A_file;
model.B_file = B_file;
model.C_file = C_file;
model.source = struct('model_dir',model_dir,'files',d.files,'shot',d.shot,'time',d.time, ...
    'description',d.description);

% operating point: zeros, overridden by the descriptor
op.coil_current = zeros(n_coils,1);
op.coil_current_given = false(n_coils,1);
op.voltage = zeros(n_coils,1);
op.voltage_given = false(n_coils,1);
op.y = zeros(p,1);
op.y_known = false(p,1);
if isfield(d,'operating_point') && isstruct(d.operating_point)
    o = d.operating_point;
    if isfield(o,'coil_current')
        [op.coil_current, op.coil_current_given] = op_values(o.coil_current,coil_names, ...
            op.coil_current,'operating_point.coil_current',desc_file);
    end
    if isfield(o,'voltage')
        [op.voltage, op.voltage_given] = op_values(o.voltage,coil_names, ...
            op.voltage,'operating_point.voltage',desc_file);
    end
    if isfield(o,'outputs')
        [op.y, op.y_known] = op_values(o.outputs,outputs,op.y,'operating_point.outputs',desc_file);
    end
end
model.op = op;

x0 = zeros(n,1);

end

function names = name_list(v, label, desc_file)
% cellstr row of unique names from a jsondecode'd JSON array of strings
if isempty(v)
    names = {};
elseif ischar(v) || isstring(v)
    names = {char(v)};
elseif iscellstr(v)
    names = v(:)';
else
    error('plasmaless:descriptor','%s: %s must be a list of names.', desc_file, label);
end
if numel(unique(names)) ~= numel(names)
    error('plasmaless:descriptor','%s: %s has duplicate names.', desc_file, label);
end
end

function G = input_matrix(d, inputs, coil_names, desc_file)
% G (m x n_coils) from the descriptor's input_map; without input_map the
% inputs must be the coil set and G is the permutation to coil_names order
m = numel(inputs); n_coils = numel(coil_names);
G = zeros(m,n_coils);
if ~isfield(d,'input_map')
    if ~isequal(sort(inputs),sort(coil_names))
        error('plasmaless:descriptor', ['%s: without input_map, inputs {%s} must be the ', ...
            'machine-description coils {%s}.'], desc_file, strjoin(inputs,','), strjoin(coil_names,','));
    end
    for jj = 1:n_coils
        G(strcmp(inputs,coil_names{jj}),jj) = 1;
    end
    return
end
if ~isstruct(d.input_map)
    error('plasmaless:descriptor','%s: input_map must be an object input -> {coil: coefficient}.', desc_file);
end
keys = fieldnames(d.input_map)';
for ii = 1:numel(keys)
    k = find(strcmp(inputs,keys{ii}),1);
    if isempty(k)
        error('plasmaless:descriptor','%s: input_map.%s is not an input.', desc_file, keys{ii});
    end
    entry = d.input_map.(keys{ii});
    if ~isstruct(entry)
        error('plasmaless:descriptor','%s: input_map.%s must be an object coil -> coefficient.', ...
            desc_file, keys{ii});
    end
    coils = fieldnames(entry)';
    for jj = 1:numel(coils)
        j = find(strcmp(coil_names,coils{jj}),1);
        if isempty(j)
            error('plasmaless:descriptor','%s: input_map.%s.%s is not a machine-description coil.', ...
                desc_file, keys{ii}, coils{jj});
        end
        v = entry.(coils{jj});
        if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v)
            error('plasmaless:descriptor','%s: input_map.%s.%s must be a finite number.', ...
                desc_file, keys{ii}, coils{jj});
        end
        G(k,j) = v;
    end
end
undriven = inputs(~any(G,2));
if ~isempty(undriven)
    warning('plasmaless:inputUndriven', ['%s: inputs {%s} have no pf_active coil voltage ', ...
        '(zero rows of input_map): not drivable through pf_active.'], desc_file, strjoin(undriven,','));
end
end

function conv = sign_convention(d, outputs, desc_file)
% convention block: state_sign, input_sign (+1 or -1), output_sign (output ->
% +1 or -1); everything absent defaults to +1 (files in the wrapper's frame)
p = numel(outputs);
conv.state_sign = 1;
conv.input_sign = 1;
conv.output_sign = ones(p,1);
if ~isfield(d,'convention')
    return
end
c = d.convention;
if ~isstruct(c)
    error('plasmaless:descriptor','%s: convention must be an object.', desc_file);
end
unknown = setdiff(fieldnames(c)', {'state_sign','input_sign','output_sign'});
if ~isempty(unknown)
    error('plasmaless:descriptor','%s: convention has unknown field(s) {%s}.', ...
        desc_file, strjoin(unknown,','));
end
if isfield(c,'state_sign')
    conv.state_sign = sign_value(c.state_sign,'convention.state_sign',desc_file);
end
if isfield(c,'input_sign')
    conv.input_sign = sign_value(c.input_sign,'convention.input_sign',desc_file);
end
if isfield(c,'output_sign')
    if ~isstruct(c.output_sign)
        error('plasmaless:descriptor','%s: convention.output_sign must be an object output -> +1/-1.', ...
            desc_file);
    end
    names = fieldnames(c.output_sign)';
    for ii = 1:numel(names)
        k = find(strcmp(outputs,names{ii}),1);
        if isempty(k)
            error('plasmaless:descriptor','%s: convention.output_sign.%s is not an output.', ...
                desc_file, names{ii});
        end
        conv.output_sign(k) = sign_value(c.output_sign.(names{ii}), ...
            ['convention.output_sign.' names{ii}], desc_file);
    end
end
end

function s = sign_value(v, label, desc_file)
% a convention entry: the number +1 or -1
if ~isnumeric(v) || ~isscalar(v) || ~(v == 1 || v == -1)
    error('plasmaless:descriptor','%s: %s must be +1 or -1.', desc_file, label);
end
s = double(v);
end

function M = load_matrix(filename)
if ~isfile(filename)
    error('plasmaless:file','Matrix file %s not found.', filename);
end
[~, ~, ext] = fileparts(filename);
switch lower(ext)
    case '.npy'
        M = read_npy(filename);
    case '.txt'
        M = readmatrix(filename,'FileType','text');
    otherwise
        error('plasmaless:descriptor','%s: matrix files must be .npy or .txt.', filename);
end
end

function [vals, given] = op_values(o, names, vals, label, desc_file)
% apply the key -> value entries of a descriptor object to vals (indexed by names)
if ~isstruct(o)
    error('plasmaless:descriptor','%s: %s must be an object name -> value.', desc_file, label);
end
given = false(size(vals));
keys = fieldnames(o)';
for ii = 1:numel(keys)
    k = find(strcmp(names,keys{ii}),1);
    v = o.(keys{ii});
    if isempty(k)
        error('plasmaless:descriptor','%s: %s.%s is not a known name.', desc_file, label, keys{ii});
    end
    if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v)
        error('plasmaless:descriptor','%s: %s.%s must be a finite number.', desc_file, label, keys{ii});
    end
    vals(k) = v;
    given(k) = true;
end
end
