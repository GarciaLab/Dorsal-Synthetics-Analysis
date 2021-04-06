function [fraction_active, mean_onset] = timesim_interp(time_vec,dls,varargin)

% numBins is the total number of bins to partition the Dorsal fluorescence;
% bin is which Dorsal fluorescence bin we're simulating

dt = 10/60; %s time resolution, should be the same as in time_vec
c = 1.5;
t_cycle = 8;
%time_vec = linspace(0, t_cycle, ceil(t_cycle/dt))';
kd = 1E3;
%dls =  kd*ones(length(ceil(t_cycle/dt)))';
nSims = 2E4;
nOffStates = 5;

% grab the Dorsal fluos over time 



%%
%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end


dls = dls(time_vec < t_cycle);
time_vec = time_vec(time_vec < t_cycle); 



q = .01;
tq = min(time_vec):q:max(time_vec);
try
    dq = makima(time_vec, dls, tq)';
catch
    dq = dls(1)*ones(length(tq), 1); 
end
% figure; plot(time_vec, dls,'o',tq,dq,'.')

onState = nOffStates + 1;


onsets = nan(nSims, 1);
for n = 1:nSims
    
    t = 0;
    state = 1;
    while t < t_cycle
        
        current_time_index = nearestIndex(tq, t); 
                
        holding_time = exprnd( (c*(dq(current_time_index)./kd) ./ (1 + dq(current_time_index)./kd) ).^-1);
        t = t + holding_time; 
        state = state + 1;
        
        if state >= onState
           onsets(n) = t;
           break;
        end

    end
    
end

onsets(onsets >= t_cycle) = nan;
fraction_active = numel(onsets(~isnan(onsets)))/nSims;
mean_onset = mean(onsets, 'omitnan');

%disp('debug stop')
end