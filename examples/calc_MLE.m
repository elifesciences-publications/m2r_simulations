% by helix length: 3, 4, 5, 6, >=7
temp_helix = {[], [], [], [], []};
% use mean ratio of all pairs in a helix
for i = 1:length(avg_results);
    idx = min(5, avg_results(i).len - 2);
    temp_helix{idx} = [temp_helix{idx}; mean(avg_results(i).ratio), mean(avg_results(i).bpp)];
end;
temp_all = [];
for i = 1:length(avg_results);
    temp_all = [temp_all; mean(avg_results(i).ratio), mean(avg_results(i).bpp)];
end;


bins = 0:0.05:1.25;
lbls = cell(1,length(bins));
for j = 1:length(bins);
    lbls{j} = ['(',num2str((j-1)*0.05),',',num2str(j*0.05),']  '];
end;
lbls{end} = [lbls{end}(1:(end-6)), 'Inf)  '];


temp_helix_bins = {};
for i = 1:5;
    temp = cell(3, length(bins));
    temp(1,:) = num2cell(bins);
    
    % hist bin by helix len
    for j = 1:size(temp_helix{i},1);
        x = find(temp_helix{i}(j,1) > bins);
        temp{2,x(end)} = [temp{2,x(end)}, temp_helix{i}(j,2)];
    end;
    
    % calculate these percentiles for helix freq ratio
    % for a given helix length
    % 2.5%, 5%, 50%, 95%, 97.5%
    for j = 1:length(bins);
        temp{3,j} = [ ...
            prctile(temp{2,j}, 2.5), ...
            prctile(temp{2,j}, 50), ...
            prctile(temp{2,j}, 97.5), ...
            ];
        
        bpp_bins = 0:0.05:1;
        bpp_temp = cell(1, length(bpp_bins));
        for k = 1:length(temp{2,j});
            x = find(temp{2,j}(k) > bpp_bins);
            bpp_temp{x(end)} = [bpp_temp{x(end)}, temp{2,j}(k)];
        end;
        for k = 1:length(bpp_bins);
            bpp_bins(k) = length(bpp_temp{k});
        end;
        [~, idx] = max(bpp_bins);
        temp{4,j} = mean(bpp_temp{idx});
        
        bpp_bins = 0:0.1:1;
        bpp_temp = cell(1, length(bpp_bins));
        for k = 1:length(temp{2,j});
            x = find(temp{2,j}(k) > bpp_bins);
            bpp_temp{x(end)} = [bpp_temp{x(end)}, temp{2,j}(k)];
        end;
        for k = 1:length(bpp_bins);
            bpp_bins(k) = length(bpp_temp{k});
        end;
        [~, idx] = max(bpp_bins);
        temp{4,j} = [temp{4,j}, mean(bpp_temp{idx})];
    end;
    
    % cell(1,:) -> bins
    % cell(2,:) -> bpp values
    % cell(3,:) -> percentiles
    % cell(4,:) -> MLE median, 0.05, 0.1
    temp_helix_bins{i} = temp;
end;

temp_all_bins = cell(4, length(bins));
temp_all_bins(1,:) = num2cell(bins);
for i = 1:size(temp_all,1);
    x = find(temp_all(i,1) > bins);
    temp_all_bins{2,x(end)} = [temp_all_bins{2,x(end)}, temp_all(i,2)];
end
% 2.5%, 5%, 50%, 95%, 97.5%
for j = 1:size(temp_all_bins,2);
    temp_all_bins{3,j} = [ ...
        prctile(temp_all_bins{2,j}, 2.5), ...
        prctile(temp_all_bins{2,j}, 50), ...
        prctile(temp_all_bins{2,j}, 97.5), ...
        ];
    
    bpp_bins = 0:0.05:1;
    bpp_temp = cell(1, length(bpp_bins));
    for k = 1:length(temp_all_bins{2,j});
        x = find(temp_all_bins{2,j}(k) > bpp_bins);
        bpp_temp{x(end)} = [bpp_temp{x(end)}, temp_all_bins{2,j}(k)];
    end;
    for k = 1:length(bpp_bins);
        bpp_bins(k) = length(bpp_temp{k});
    end;
    [~, idx] = max(bpp_bins);
    temp_all_bins{4,j} = mean(bpp_temp{idx});
    
    bpp_bins = 0:0.1:1;
    bpp_temp = cell(1, length(bpp_bins));
    for k = 1:length(temp_all_bins{2,j});
        x = find(temp_all_bins{2,j}(k) > bpp_bins);
        bpp_temp{x(end)} = [bpp_temp{x(end)}, temp_all_bins{2,j}(k)];
    end;
    for k = 1:length(bpp_bins);
        bpp_bins(k) = length(bpp_temp{k});
    end;
    [~, idx] = max(bpp_bins);
    temp_all_bins{4,j} = [temp_all_bins{4,j}, mean(bpp_temp{idx})];
end;

