function dict_to_smr(data, fname, spk_freq, lfp_freq)
fhand = CEDS64Create( fname, 32, 1 );
CEDS64TimeBase( fhand, 1./spk_freq );

names = fieldnames(data);

for i = 1:length(names)
    curr_name = char(names(i));
    if startsWith(curr_name, 'lfp', 'IgnoreCase', true)
        CEDS64SetWaveChan( fhand, i, 32, 9, lfp_freq );  
    else
        CEDS64SetWaveChan( fhand, i, 1, 9, spk_freq);
    end
    CEDS64WriteWave( fhand, i, data.(curr_name), 0);
    CEDS64ChanTitle( fhand, i, curr_name );
    CEDS64ChanUnits( fhand, i, 'V' );
end

CEDS64CloseAll();
end