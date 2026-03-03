function bytes = get_approx_free_bytes()
    bytes = [];
    try
        fid = fopen('/proc/meminfo','r');
        if fid > 0
            s = textscan(fid, '%s %s %s', 'Delimiter', ':');
            fclose(fid);
            keys = s{1};
            vals = s{2};
            ii = find(strcmp(keys,'MemAvailable') | strcmp(keys,'MemFree') | strcmp(keys,'MemTotal'), 1);
            if ~isempty(ii)
                v = vals{ii};
                % v like '12345678 kB'
                tok = textscan(v, '%f %s');
                kb = tok{1};
                bytes = single(kb * 1024);
            end
        end
    catch
        bytes = [];
    end
end
