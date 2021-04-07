function full_theta = getAllThetas(full_names, results)

%Get the mean values for the optimized parameters (from results.mean) and the default 
%values for the unoptimized parameters (from results.theta, the values from
%the very last mcmc step)

%full_names = ["c", "kd" , "nentrystates", "moffstates", "pentry", "pexit", "tcycle"];
chain_names = string(results.names);
batched_indices = find(contains(chain_names, '['));
batched_names = chain_names(contains(chain_names, '['));
n_batches = length(batched_names);
batched_name = extractBefore(batched_names(2), '[');
batched_full_index = find(full_names == batched_name);


full_theta = cell(1, n_batches);
for j = 1:n_batches
    full_theta{j} = nan(1, length(full_names));
    full_theta{j}(batched_full_index) = results.mean(batched_indices(j));
end

for k = 1:length(full_names)
    if k ~= batched_full_index
        chain_ind = find(full_names(k) == chain_names);
        if ~isempty(chain_ind)
            for j = 1:length(full_theta)
                full_theta{j}(k) = results.mean(chain_ind);
            end
        else
            for j = 1:length(full_theta)
                %the RHS of this will have the wrong index if 
                %the batched parameter is not the second parameter
                %(kd). 
                full_theta{j}(k) = results.theta(k + n_batches - 1); 
            end
        end
    end
end