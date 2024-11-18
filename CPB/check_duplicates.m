function [time_new signal_new] = check_duplicates(time, signal)

NDP = numel(signal);
p = 1;
index_remove = 0;
for n = 2:NDP
    
    if signal(n)==signal(n-1)
        index_remove(p) = n;
        p = p + 1;
    else
        continue
    end
end

time_new = time;
signal_new = signal;

time_new(index_remove) = [];
signal_new(index_remove) = [];
