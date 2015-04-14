function [ panel ] = var_data_gen( len )
%VAR_DATA_GEN Generate VAR data
%

dim = 2;
lb_pos = [1 24 25]; % look-back positions
para = {[0.8 0; 0 0.8], [1 0; 0, 1], [-0.8 0; 0 -0.8]}; % AR matrices
cov = [100 30; 30 100];
mu = [100; 150];
%mu = [0; 0];


init_y = [80 120; ...
          100 150; 100 150; 120 180; 120 180; 140 210; 140 210; ...
          160 240; 160 240; 140 210; 140 210; 120 180; 120 180; ...
          100 150; 100 150; 80 120; 80 120; 60 90; 60 90; ...
          40 60; 40 60; 60 90; 60 90; 80 120; 80 120 ]'; 
      % data for init, must have dim * lb dimensions

lb = max(lb_pos);
panel = zeros(dim, len + lb);
panel(:, 1:lb) = init_y;
for ii = 1:len
    panel(:, lb+ii) = mvnrnd(zeros(2, 1), cov, 1)' + mu;
    for ipos = 1:length(lb_pos)
        panel(:, lb+ii) = panel(:, lb+ii) + ...
            para{ipos} * (panel(:, lb+ii-lb_pos(ipos)) - mu);
    end
end

