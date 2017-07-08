% function dict_to_smr(data, fname, spk_freq, lfp_freq)
% disp(fieldnames(data));
% 
% fhand2 = CEDS64Create( fname, 32, 1 );
% CEDS64TimeBase( fhand2, 1./spk_freq );
% 
% CEDS64SetWaveChan( fhand2, 1, 1, 9, spk_freq );
% [ i64TNext ] = CEDS64WriteWave( fhand2, 1, data.Spk_10033, 0);
% disp(i64TNext);
% 
% CEDS64SetWaveChan( fhand2, 2, 1, 9, spk_freq );
% [ i64TNext ] = CEDS64WriteWave( fhand2, 2, data.Spk_10034, 0);
% disp(i64TNext);
% 
% CEDS64SetWaveChan( fhand2, 3, 1, 9, spk_freq );
% [ i64TNext ] = CEDS64WriteWave( fhand2, 3, data.Spk_10035, 0);
% disp(i64TNext);
% 
% CEDS64SetWaveChan( fhand2, 4, 32, 9, lfp_freq );
% [ i64TNext ] = CEDS64WriteWave( fhand2, 4, data.LFP_10006, 0);
% disp(i64TNext);
% 
% CEDS64SetWaveChan( fhand2, 5, 32, 9, lfp_freq );
% [ i64TNext ] = CEDS64WriteWave( fhand2, 5, data.LFP_10007, 0);
% disp(i64TNext);
% 
% CEDS64SetWaveChan( fhand2, 6, 32, 9, lfp_freq );
% [ i64TNext ] = CEDS64WriteWave( fhand2, 6, data.LFP_10008, 0);
% disp(i64TNext);
% 
% CEDS64CloseAll();
% 
% end

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