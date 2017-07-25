MLE = load('add_20170716_RFAM_MLE_95.mat', 'temp_helix_bins');
load('add_20160811_a140T_score_cleanup_more.mat');

% weights = ones(1,128);
% weights = [ones(1,7)*0.5, ones(1,10), ones(1,4)*0.5, ones(1,7), 0.5, ones(1,94), 0.5, ones(1,4)];
weights = [ones(1,74), 0,0,0, ones(1,51)];
% weights = [ones(1,7)*0.5, ones(1,10), ones(1,4)*0.5, ones(1,7), 0.5, ones(1,45), 0,0,0, ones(1,46), 0.5, ones(1,4)];

for i = 1:length(quartet_MiP1_exp);
    quartet_MiP1_exp{i}(quartet_MiP1_exp{i}>1)=1;
    dist_MiP1_exp(i,:) = get_quartet_dist(quartet_MiP1_exp{i}, weights);
end
for i = 1:length(quartet_MiP2_exp);
    quartet_MiP2_exp{i}(quartet_MiP2_exp{i}>1)=1;
    dist_MiP2_exp(i,:) = get_quartet_dist(quartet_MiP2_exp{i}, weights);
end
for i = 1:length(quartet_MiP4a_exp);
    quartet_MiP4a_exp{i}(quartet_MiP4a_exp{i}>1)=1;
    dist_MiP4a_exp(i,:) = get_quartet_dist(quartet_MiP4a_exp{i}, weights);
end
for i = 1:length(quartet_MiP4a_more_exp);
    quartet_MiP4a_more_exp{i}(quartet_MiP4a_more_exp{i}>1)=1;
    dist_MiP4a_more_exp(i,:) = get_quartet_dist(quartet_MiP4a_more_exp{i}, weights);
end
for i = 1:length(quartet_MiP4b_exp);
    quartet_MiP4b_exp{i}(quartet_MiP4b_exp{i}>1)=1;
    dist_MiP4b_exp(i,:) = get_quartet_dist(quartet_MiP4b_exp{i}, weights);
end
for i = 1:length(quartet_MtP1b_exp);
    quartet_MtP1b_exp{i}(quartet_MtP1b_exp{i}>1)=1;
    dist_MtP1b_exp(i,:) = get_quartet_dist(quartet_MtP1b_exp{i}, weights);
end


% lock p1
MiP1 = dist_MiP1_exp(:,1) ./ max(dist_MiP1_exp(:,2),dist_MiP1_exp(:,3));
get_bpp_estimate(MLE.temp_helix_bins, MiP1(1:6)); % p2
get_bpp_estimate(MLE.temp_helix_bins, MiP1(7:10)); % p4b
% lock p2
MiP2 = dist_MiP2_exp(:,1) ./ max(dist_MiP2_exp(:,2),dist_MiP2_exp(:,3));
get_bpp_estimate(MLE.temp_helix_bins, MiP2(1:7)); % p1
get_bpp_estimate(MLE.temp_helix_bins, MiP2(10:13)); % p4a
% lock p4a
MiP4A = dist_MiP4a_exp(:,1) ./ max(dist_MiP4a_exp(:,2),dist_MiP4a_exp(:,3));
get_bpp_estimate(MLE.temp_helix_bins, MiP4A(1:5)); % p1b
get_bpp_estimate(MLE.temp_helix_bins, MiP4A(6:7)); % p2b
MiP4A_more = dist_MiP4a_more_exp(:,1) ./ max(dist_MiP4a_more_exp(:,2),dist_MiP4a_more_exp(:,3));
get_bpp_estimate(MLE.temp_helix_bins, MiP4A_more); % p2
% lock p4b
MiP4B = dist_MiP4b_exp(:,1) ./ max(dist_MiP4b_exp(:,2),dist_MiP4b_exp(:,3));
get_bpp_estimate(MLE.temp_helix_bins, MiP4B(18:22)); % p1b
get_bpp_estimate(MLE.temp_helix_bins, MiP4B(23:24)); % p2b
get_bpp_estimate(MLE.temp_helix_bins, MiP4B(1:7)); % p1
get_bpp_estimate(MLE.temp_helix_bins, MiP4B(8:13)); % p2
get_bpp_estimate(MLE.temp_helix_bins, MiP4B(14:17)); % p4a
% lock p1b
MtP1B = dist_MtP1b_exp(:,1) ./ max(dist_MtP1b_exp(:,2),dist_MtP1b_exp(:,3));
get_bpp_estimate(MLE.temp_helix_bins, MtP1B(9:10)); % p2b
get_bpp_estimate(MLE.temp_helix_bins, MtP1B(1:4)); % p4a
get_bpp_estimate(MLE.temp_helix_bins, MtP1B(5:8)); % p4b
get_bpp_estimate(MLE.temp_helix_bins, MtP1B(11:16)); % p3

% mut p2
MtP2 = dist_MtP2_exp(:,1) ./ max(dist_MtP2_exp(:,2),dist_MtP2_exp(:,3));
get_bpp_estimate(MLE.temp_helix_bins, MtP2(1:7)); % p1
get_bpp_estimate(MLE.temp_helix_bins, MtP2(10:13)); % p4a


% (-) add M2R
rdat_M2R.area_norm_nolig(rdat_M2R.area_norm_nolig>1)=1;
for i = 1:95;
    quartet_M2R_exp{i} = rdat_M2R.area_norm_nolig(27:154, rdat_M2R.all{i}.index)';
    quartet_M2R_exp{i}(quartet_M2R_exp{i}>1)=1;
    dist_M2R_exp(i,:) = get_quartet_dist(quartet_M2R_exp{i}, weights);
end
M2R_nolig = dist_M2R_exp(:,1) ./ max(dist_M2R_exp(:,2),dist_M2R_exp(:,3));

get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(1:7)); % p1
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(8:13)); % p2
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(14:19)); % p3
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(20:25)); % p5
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(26:30)); % p1b
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(31:32)); % p2b
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(33:36)); % p4a
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(37:40)); % p4b
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(41:42)); % p4c
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(43:48)); % p6
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(49:52)); % p8
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(53:57)); % p9
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(58:60)); % p10
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(61:69)); % p11
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(70:72)); % p12
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(73:76)); % p13
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(77:79)); % p14
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(80:86)); % p15
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(87:91)); % p16
get_bpp_estimate(MLE.temp_helix_bins, M2R_nolig(92:95)); % p17

% (+) add M2R
rdat_M2R.area_norm_lig(rdat_M2R.area_norm_lig>1)=1;
for i = 1:95;
    quartet_M2R_exp_lig{i} = rdat_M2R.area_norm_lig(27:154, rdat_M2R.all{i}.index)';
    quartet_M2R_exp_lig{i}(quartet_M2R_exp_lig{i}>1)=1;
    dist_M2R_exp_lig(i,:) = get_quartet_dist(quartet_M2R_exp_lig{i}, weights);
end
M2R_lig = dist_M2R_exp_lig(:,1) ./ max(dist_M2R_exp_lig(:,2),dist_M2R_exp_lig(:,3));

get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(1:7)); % p1
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(8:13)); % p2
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(14:19)); % p3
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(20:25)); % p5
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(26:30)); % p1b
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(31:32)); % p2b
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(33:36)); % p4a
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(37:40)); % p4b
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(41:42)); % p4c
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(43:48)); % p6
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(49:52)); % p8
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(53:57)); % p9
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(58:60)); % p10
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(61:69)); % p11
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(70:72)); % p12
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(73:76)); % p13
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(77:79)); % p14
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(80:86)); % p15
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(87:91)); % p16
get_bpp_estimate(MLE.temp_helix_bins, M2R_lig(92:95)); % p17


%% Ann's result anyway
rdat_supp = load('add_20161014_3Dlock_quartet');

pairs_P4a_supp = zeros(2,4);
temp = rdat_supp.keys_P4a(4:3:end);
for i = 1:length(temp);
    pos = strsplit(temp{i}, ';');
    pairs_P4a_supp(1,i) = str2num(pos{1}(2:end-1));
    pairs_P4a_supp(2,i) = str2num(pos{2}(2:end-1));
end;
pairs_P2_supp = zeros(2,4);
temp = rdat_supp.keys_P2(4:3:end);
for i = 1:length(temp);
    pos = strsplit(temp{i}, ';');
    pairs_P2_supp(1,i) = str2num(pos{1}(2:end-1));
    pairs_P2_supp(2,i) = str2num(pos{2}(2:end-1));
end;

rdat_supp.d_MiP4a_nolig(rdat_supp.d_MiP4a_nolig>1)=1;
rdat_supp.d_MtP4a_nolig(rdat_supp.d_MtP4a_nolig>1)=1;
rdat_supp.d_MiP2_nolig(rdat_supp.d_MiP2_nolig>1)=1;
rdat_supp.d_MtP2_nolig(rdat_supp.d_MtP2_nolig>1)=1;

for i =1:4
    quartet_MiP4a_supp_exp{i} = rdat_supp.d_MiP4a_nolig(27:154,[1, ((i-1)*3+2):i*3+1])';
    dist_MiP4a_supp_exp(i,:) = get_quartet_dist(quartet_MiP4a_supp_exp{i}, weights);
end
for i =1:4
    quartet_MtP4a_supp_exp{i} = rdat_supp.d_MtP4a_nolig(27:154,[1, ((i-1)*3+2):i*3+1])';
    dist_MtP4a_supp_exp(i,:) = get_quartet_dist(quartet_MtP4a_supp_exp{i}, weights);
end
for i =1:4
    quartet_MiP2_supp_exp{i} = rdat_supp.d_MiP2_nolig(27:154,[1, ((i-1)*3+2):i*3+1])';
    dist_MiP2_supp_exp(i,:) = get_quartet_dist(quartet_MiP2_supp_exp{i}, weights);
end
for i =1:4
    quartet_MtP2_supp_exp{i} = rdat_supp.d_MtP2_nolig(27:154,[1, ((i-1)*3+2):i*3+1])';
    dist_MtP2_supp_exp(i,:) = get_quartet_dist(quartet_MtP2_supp_exp{i}, weights);
end

% lock p4a
MiP4A_supp = dist_MiP4a_supp_exp(:,1) ./ max(dist_MiP4a_supp_exp(:,2),dist_MiP4a_supp_exp(:,3));
get_bpp_estimate(MLE.temp_helix_bins, MiP4A_supp); % p4b
% lock p2
MiP2_supp = dist_MiP2_supp_exp(:,1) ./ max(dist_MiP2_supp_exp(:,2),dist_MiP2_supp_exp(:,3));
get_bpp_estimate(MLE.temp_helix_bins, MiP2_supp); % p4b
% mut p2
MtP2_supp = dist_MtP2_supp_exp(:,1) ./ max(dist_MtP2_supp_exp(:,2),dist_MtP2_supp_exp(:,3));
get_bpp_estimate(MLE.temp_helix_bins, MtP2_supp); % p4b



%% 2-bp M2R
for i = 1:44;
    quartet_3Dst_exp{i} = rdat_3Dst.area_norm_nolig(27:154,[1, ((i-1)*3+2):i*3+1])';
    quartet_3Dst_exp{i}(quartet_3Dst_exp{i}>1)=1;
    dist_3Dst_exp(i,:) = get_quartet_dist(quartet_3Dst_exp{i}, weights);
end
for i = 1:44;
    quartet_3Dst_exp_lig{i} = rdat_3Dst.area_norm_lig(27:154,[1, ((i-1)*3+2):i*3+1])';
    quartet_3Dst_exp_lig{i}(quartet_3Dst_exp_lig{i}>1)=1;
    dist_3Dst_exp_lig(i,:) = get_quartet_dist(quartet_3Dst_exp_lig{i}, weights);
end
D3st_nolig = dist_3Dst_exp(:,1) ./ max(dist_3Dst_exp(:,2),dist_3Dst_exp(:,3));
D3st_lig = dist_3Dst_exp_lig(:,1) ./ max(dist_3Dst_exp_lig(:,2),dist_3Dst_exp_lig(:,3));

% nolig
get_bpp_estimate(MLE.temp_helix_bins, D3st_nolig(1:6), 1); % p1
get_bpp_estimate(MLE.temp_helix_bins, D3st_nolig(7:11), 1); % p2
get_bpp_estimate(MLE.temp_helix_bins, D3st_nolig(12:16), 1); % p3
get_bpp_estimate(MLE.temp_helix_bins, D3st_nolig(17:21), 1); % p5
get_bpp_estimate(MLE.temp_helix_bins, D3st_nolig(22:25), 1); % p1b
get_bpp_estimate(MLE.temp_helix_bins, D3st_nolig(26), 1); % p2b
get_bpp_estimate(MLE.temp_helix_bins, D3st_nolig(29:31), 1); % p4a
get_bpp_estimate(MLE.temp_helix_bins, D3st_nolig(32:34), 1); % p4b
get_bpp_estimate(MLE.temp_helix_bins, D3st_nolig(35), 1); % p4c
get_bpp_estimate(MLE.temp_helix_bins, D3st_nolig(36:40), 1); % p6
get_bpp_estimate(MLE.temp_helix_bins, D3st_nolig(41:43), 1); % p13
% lig
get_bpp_estimate(MLE.temp_helix_bins, D3st_lig(1:6), 1); % p1
get_bpp_estimate(MLE.temp_helix_bins, D3st_lig(7:11), 1); % p2
get_bpp_estimate(MLE.temp_helix_bins, D3st_lig(12:16), 1); % p3
get_bpp_estimate(MLE.temp_helix_bins, D3st_lig(17:21), 1); % p5
get_bpp_estimate(MLE.temp_helix_bins, D3st_lig(22:25), 1); % p1b
get_bpp_estimate(MLE.temp_helix_bins, D3st_lig(26), 1); % p2b
get_bpp_estimate(MLE.temp_helix_bins, D3st_lig(29:31), 1); % p4a
get_bpp_estimate(MLE.temp_helix_bins, D3st_lig(32:34), 1); % p4b
get_bpp_estimate(MLE.temp_helix_bins, D3st_lig(35), 1); % p4c
get_bpp_estimate(MLE.temp_helix_bins, D3st_lig(36:40), 1); % p6
get_bpp_estimate(MLE.temp_helix_bins, D3st_lig(41:43), 1); % p13

