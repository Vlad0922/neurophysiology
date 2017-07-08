function dict_to_smr(data, fname, spk_freq, lfp_freq)
fhand2 = CEDS64Create( fname, 32, 1 );
CEDS64TimeBase( fhand2, 1./spk_freq );

names = fieldnames(data);

for i = 1:length(names)
    curr_name = char(names(i));
    if startsWith(curr_name, 'lfp', 'IgnoreCase', true)
        CEDS64SetWaveChan( fhand2, i, 32, 9, lfp_freq );  
    else
        CEDS64SetWaveChan( fhand2, i, 1, 9, spk_freq);
    end
    CEDS64WriteWave( fhand2, i, data.(curr_name), 0);
    CEDS64ChanTitle( fhand2, i, curr_name );
end

CEDS64CloseAll();
end