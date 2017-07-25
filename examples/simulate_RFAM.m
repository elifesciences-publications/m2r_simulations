% Store results
struct_is = [];

% To start over, use 1:x [maybe 20 at a time]
% Run in batches since it takes a very long time
for n = 334:350;
    % skip RFAM sequences that has "N" in it
    if ~isempty(strrep(strrep(strrep(strrep(rf_seq{n}, 'A', ''), 'C', ''), 'G', ''), 'U', '')); continue; end;
    struct_is = [struct_is, simulate_rfam_entry(rf_seq{n}, rf_name{n})];
end;


%%% by helix length
% min helix length, i.e. only consider helices lenght >= 3
N_SAMPLE = 2; 
% min bpp cutoff, i.e. helices with bpp in WT < 0.01 are discarded
N_CUTOFF = 0.01;

global_idx = 0;
avg_results = [];
for i = 1:length(struct_is);
    for j = 1:length(struct_is(i).helix_len);
        % get indices of bpp-valid helices
        temp_idx = [];
        for k = 1:struct_is(i).helix_len(j);
            if bpp_all(global_idx + k) > N_CUTOFF;
                temp_idx = [temp_idx, k];
            end;
        end;
        % increment global index, skip short helices
        if length(temp_idx) <= N_SAMPLE;
            global_idx = global_idx + struct_is(i).helix_len(j);
            continue;
        end;
        
        temp_dist = [];
        temp_bpp = [];
        temp_ratio = [];
        for k = 1:length(temp_idx);
            temp_dist = [temp_dist; dist_all(global_idx + temp_idx(k),:)];
            temp_bpp = [temp_bpp, bpp_all(global_idx + temp_idx(k))];
            % the helix frequency "ratio"
            temp_ratio = [temp_ratio, dist_all(global_idx + temp_idx(k),1) / max(dist_all(global_idx + temp_idx(k),2), dist_all(global_idx + temp_idx(k),3))];
        end;
        
        temp = [];
        temp.bpp = temp_bpp;
        temp.dist = temp_dist;
        temp.ratio = temp_ratio;
        temp.len = length(temp_idx);
        % flatten all helices of all RFAM families into one giant list
        avg_results = [avg_results, temp];
        
        global_idx = global_idx + struct_is(i).helix_len(j);
    end;
end;
