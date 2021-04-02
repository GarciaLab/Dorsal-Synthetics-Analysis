function [fraction_active, mean_onset] = timesim(varargin)

dt = 10/60; %s time resolution
c = 1;
t_cycle = 8;
kd = 1E3;
dls =  kd*ones(length(time_vec), 1);
nSims = 20;
nOffStates = 5;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

time_vec = linspace(0, t_cycle, ceil(t_cycle/dt))';

onState = nOffStates + 1;


onsets = nan(nSims, 1);

for n = 1:nSims
    
    state = ones(length(time_vec), 1);
    
    tauoffs = exprnd( (c*(dls./kd) ./ (1 + dls./kd) ).^-1, [length(dls), 1]);
    
    for t = 1:length(time_vec)
        
        if tauoffs(t) < dt
            state(t:end) = state(t) + 1;
        end
        
        if state(t) == onState
            onsets(n) = time_vec(t);
            break;
        end
        
    end
    
    
end

fraction_active = numel(onsets(~isnan(onsets)))/nSims;
mean_onset = mean(onsets, 'omitnan');

disp('debug stop')
end