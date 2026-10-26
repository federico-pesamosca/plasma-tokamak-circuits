function key = plasmaless_supply_key(supply)
% short power-supply key (e.g. 'CS3U', 'VS1', 'VS3') of a pf_active.supply{j}
% element, from its name: a long name 'Long description (KEY)' gives 'KEY'
% (same rule as plasmaless_coil_key); otherwise the part before the first
% '-' of the machine-description names ('CS3U-1-2' -> 'CS3U',
% 'VS1-1-2-3-4-5-6' -> 'VS1', 'VS3-1' -> 'VS3'); a plain key ('CS3U') is
% returned as is.
% Shared by plasmaless_supply_names (model.supply_names) and
% plasmaless_input_voltages (matching of the input supplies).

key = strtrim(supply.name);
tok = regexp(key,'\(([^()]*)\)','tokens','once');
if ~isempty(tok)
    key = tok{1};
else
    key = regexprep(key,'-.*$','');
end
key = strtrim(key);
if isempty(key)
    error('plasmaless:supplyName','pf_active supply without name.');
end

end
