%% Concentric-shell measurements on histology (nearest-structure method)
%{
Applies the OCT shell method (measure_parenchyma_shells.m) to the
histology sections and writes one long-format CSV for the Bayesian models.

Output: one row per stain x section x object x 40 µm shell, with columns
  stain, subject, severity, region, section, object_id, obj, class,
  d_inner_um, d_outer_um, n_px, raw, hmatch, z,
  obj_area_um2, obj_diam_um, enough_pixels

Requires on path: import_stain.m, measure_histo_shells.m
%}

clc; close all; clearvars -except caa swp_op
addpath(fileparts(mfilename('fullpath')));

%% Settings
hdir = '/projectnb/npbssmic/ns/CAA/histology/';

cfg.pix_um                 = 1.3067;     % µm per pixel
cfg.edges_um               = 0:40:520;   % same shells as the OCT analysis
cfg.n_min_obj_px           = 0;          % manual annotations: keep all
cfg.n_min_bin_px           = 100;        % flag object x shell with < 100 px
cfg.epvs_buffer_px         = 1;          % vessels touching an EPVS are not "normal"
cfg.hist_match             = true;       % match to first section of each stain
cfg.zscore_parenchyma_only = true;       % z-score from parenchyma, not EPVS/vessel pixels
cfg.save_qc                = true;

% Stains and filename suffixes (as in the ring script)
stains = struct( ...
  'name',         {'LHE',                   'CD68',               'GFAP'}, ...
  'dir',          {'LHE/',                  'CD68/',              'GFAP/'}, ...
  'stain_suffix', {'_M_BackgroundMask.tif', '_shrunk_deconv.tif', '_shrunk_deconv.tif'}, ...
  'epvs_suffix',  {'_EPVS_Mask.tif',        '_shrunk_EPVS.tif',   '_shrunk_EPVS.tif'}, ...
  'ves_suffix',   {'_VESSEL_Mask.tif',      '_shrunk_VESSEL.tif', '_shrunk_VESSEL.tif'}, ...
  'mask_suffix',  {'_mask.tif',             '_shrunk_mask.tif',   '_shrunk_mask.tif'});

% Subject -> severity (same as the OCT analysis)
severity_map = containers.Map( ...
    {'caa6',    'caa26',   'caa17',    'caa22',  'caa25'}, ...
    {'control', 'control', 'moderate', 'severe', 'severe'});

% Region is read from the section name: _1 = frontal, _7 = occipital (right
% after the subject number), else "front"/"occip" in the name. For sections
% whose name does not contain the region, add it here, e.g.
%   region_override('CAA26_1') = 'front';
region_override = containers.Map('KeyType', 'char', 'ValueType', 'char');

%% Measure each stain
allT = {};
for s = 1:numel(stains)
    st = stains(s);
    fprintf('\n===== %s =====\n', st.name);
    H = import_stain(fullfile(hdir, st.dir), st.stain_suffix, st.epvs_suffix, ...
                     st.ves_suffix, st.mask_suffix);

    cfg.qc_dir = fullfile(hdir, 'shell_qc', st.name);
    T = measure_histo_shells(H, st.name, cfg);
    if isempty(T), warning('%s: no rows produced', st.name); continue; end

    %%% Metadata from the section name
    secs = unique(T.section);
    subj = strings(numel(secs), 1); reg = subj; sev = subj;
    for k = 1:numel(secs)
        [subj(k), reg(k)] = parse_section_name(secs(k), region_override);
        if isKey(severity_map, char(subj(k)))
            sev(k) = severity_map(char(subj(k)));
        else
            sev(k) = "unknown";
        end
    end
    [~, loc] = ismember(T.section, secs);
    T.stain    = repmat(string(st.name), height(T), 1);
    T.subject  = subj(loc);
    T.region   = reg(loc);
    T.severity = sev(loc);
    T.object_id = T.stain + "_" + T.section + "_" + T.class + "_" + string(T.obj);

    fprintf('\n%s sections:\n', st.name);
    disp(table(secs, subj, reg, sev, 'VariableNames', {'section','subject','region','severity'}));
    allT{end+1} = T; %#ok<SAGROW>
end

%% Write CSV
H_all = vertcat(allT{:});
H_all = movevars(H_all, {'stain','subject','severity','region','section','object_id'}, 'Before', 1);

n_unknown = nnz(H_all.region == "unknown" | H_all.severity == "unknown");
if n_unknown > 0
    warning(['%d rows have unknown region or severity. Check the section ' ...
             'tables above and add entries to region_override / severity_map.'], n_unknown);
end

dt   = string(datetime('now', 'Format', 'yyyy-MM-dd'));
fout = fullfile(hdir, "histo_shells_" + dt + ".csv");
writetable(H_all, fout);
fprintf('\nWrote %d rows to %s\n', height(H_all), fout);

%% ------------------------------------------------------------------------
function [subj, region] = parse_section_name(name, region_override)
% Subject from "CAA26", "CAA_26", "caa26" ...
% Region from the code right after the subject number:
%   _1_ or _1 at the end -> frontal;  _7_ or _7 at the end -> occipital
%   e.g. CAA26_1_GFAP_10x -> caa26 / front,  CAA_22_7 -> caa22 / occip
% Falls back to "front"/"occip" in the name, then to region_override.
n = lower(char(name));
tok = regexp(n, 'caa_?(\d+)', 'tokens', 'once');
if isempty(tok), subj = "unknown"; else, subj = "caa" + string(tok{1}); end

code = regexp(n, 'caa_?\d+_([17])(?=_|\.|$)', 'tokens', 'once');
if ~isempty(code) && strcmp(code{1}, '1')
    region = "front";
elseif ~isempty(code) && strcmp(code{1}, '7')
    region = "occip";
elseif contains(n, 'front')
    region = "front";
elseif contains(n, 'occip')
    region = "occip";
else
    region = "unknown";
    keys_ = keys(region_override);
    [~, order] = sort(cellfun(@numel, keys_), 'descend');   % longest prefix wins
    for k = 1:numel(keys_)
        key = keys_{order(k)};
        if startsWith(n, lower(key))
            region = string(region_override(key));
            break
        end
    end
end
end