function sec_since_epoch = convertToUnixEpoch(time, tz)

unix_epoch = 719529;

if tz
    % Convert from EST to UTC
    time = time + (-tz/24);
end

sec_since_epoch = 86400 * (time - unix_epoch);