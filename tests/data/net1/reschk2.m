function [ ] = reschk2( dat)
%RESCHK display 7 time plots for estep results

N_CUST = 7;
N_STEP = 20;
color = {'y', 'g', 'b', 'r', 'k'};

pop_size = size(dat, 1); % number of rows
step_pop_size  = floor(pop_size / N_STEP);

axs = cell(7, 1);

fig_1 = figure('Name', ['MCMC E-Step: trend of the population']);
for iax = 1:N_CUST
    axs{iax} = subplot(2, 4, iax);
    
    set(axs{iax}, 'fontsize', 7);
    set(axs{iax},'XTick',[]);
    set(axs{iax},'ylim',[0 340], 'ytick', 0:100:300);
    
    
    temp_dat = nan(pop_size, N_STEP);
    for itime = 1:N_STEP
        temp_dat(1:itime*step_pop_size, itime) = ...
            dat(1:itime*step_pop_size, iax);
        
    boxplot(axs{iax}, temp_dat, 'outliersize', 2, 'jitter', 0.5);
    set(axs{iax},'XTickLabel',{' '});
    %set(axs{iax},'ylim',[0 340], 'ytick', 0:100:300, 'fontsize', 8);
    line(1:N_STEP, mean(temp_dat, 1)', 'color', 'g', 'parent', axs{iax});
    
    end
end

fig_2 = figure('Name', '1D and 2D histograms');
for ii = 1:N_CUST
    for jj = 1:N_CUST
        cur_ax = subplot(N_CUST, N_CUST, (ii-1)*N_CUST+jj);
        set(cur_ax, 'fontsize', 6);
        if ii == jj
            hist(dat(:, ii));
        else
            ndhist(dat(:, ii), dat(:, jj));
        end
    end
end
%tx = 0:24:168;
%set(axs{N_CUST}, 'xtick', tx, 'xticklabel', cellstr(num2str(tx')));
    
    


end

