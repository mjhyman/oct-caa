function T = add_percent_change(T, parench_median, mus_upper, ret_upper)
% ADD_PERCENT_CHANGE  Add percentage change from the subject/region
% white-matter median, plus raw-value bound flags, to a shells table.
%
%   pct = (value - median) / median * 100
%
% INPUTS:
%   T              (table):  output of measure_parenchyma_shells.m; needs
%                            columns subject, region, mus, ret
%   parench_median (struct): parench_median.(subject).(region).med.mus/.ret
%   mus_upper      (double): raw mus bound (old pipeline: 25)
%   ret_upper      (double): raw retardance bound (old pipeline: 45)
% OUTPUTS (columns added or overwritten):
%   mus_med, ret_med : median used for each row
%   mus_pct, ret_pct : percentage change from that median
%   mus_ok,  ret_ok  : raw value below the bound (filter on these in R)

T.subject = string(T.subject);
T.region  = string(T.region);
n = height(T);
T.mus_med = NaN(n, 1);
T.ret_med = NaN(n, 1);

keys = unique(T(:, {'subject', 'region'}), 'rows');
for i = 1:height(keys)
    sub  = char(keys.subject(i));
    loc  = char(keys.region(i));
    rows = T.subject == keys.subject(i) & T.region == keys.region(i);

    if ~isfield(parench_median, sub) || ~isfield(parench_median.(sub), loc)
        warning('No median for %s %s: percentage change left as NaN', sub, loc);
        continue
    end
    m = parench_median.(sub).(loc).med;
    assert(isscalar(m.mus) && isscalar(m.ret), ...
           'Median for %s %s is not scalar', sub, loc);

    T.mus_med(rows) = m.mus;
    T.ret_med(rows) = m.ret;
    fprintf('  %s %s: median mus = %.3f, median ret = %.3f (%d rows)\n', ...
            sub, loc, m.mus, m.ret, nnz(rows));
end

T.mus_pct = (T.mus - T.mus_med) ./ T.mus_med * 100;
T.ret_pct = (T.ret - T.ret_med) ./ T.ret_med * 100;

% Bounds apply to RAW values (the old script applied them after conversion
% to percent, which cut off percentage changes above +25% / +45%)
T.mus_ok = T.mus < mus_upper;
T.ret_ok = T.ret < ret_upper;
end