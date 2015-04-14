function [ ] = reschk2( dat)
%RESCHK display 7 time plots for estep results

N_CUST = 7;
N_STEP = 20;
color = {'y', 'g', 'b', 'r', 'k'};


fig_2 = figure('Name', '1D and 2D histograms');
for ii = 1:N_CUST
    for jj = 1:N_CUST
        cur_ax = subplot(N_CUST, N_CUST, (ii-1)*N_CUST+jj);
        set(cur_ax, 'fontsize', 6);
        if ii == jj
            hist(dat(:, ii));
        else
            ndhist(dat(:, ii), dat(:, jj));
            sf = findobj(cur_ax, 'type', 'surface');
            set(sf, 'zdata', get(sf, 'cdata'));
        end
    end
end
%tx = 0:24:168;
%set(axs{N_CUST}, 'xtick', tx, 'xticklabel', cellstr(num2str(tx')));
    
    


end

