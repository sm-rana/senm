function [ ] = reschk( dat, win)
%RESCHK display 7 time plots for estep results

N_CUST = 7;
N_TIME = 168;
N_STEP = 4;
color = {'y', 'g', 'b', 'r', 'k'};

pop_size = size(dat, 1); % number of rows

fig = figure('Name', ['MCMC E-Step: ' num2str(pop_size), ' samples']);

axs = cell(7, 1);

for iax = 1:N_CUST
    axs{iax} = subplot(2, 4, iax);
    part_dat = dat(:, end-win*N_CUST+iax:N_CUST:end);
    boxplot(axs{iax}, part_dat, 'outliersize', 1);
    set(axs{iax},'XTickLabel',{' '});
    set(axs{iax},'ylim',[0 340], 'ytick', 0:100:300, 'fontsize', 8);
    line(1:win, mean(part_dat, 1)', 'color', 'g', 'parent', axs{iax});
    
    
    %set(axs{iax},'xtickmode','auto','xticklabelmode','auto');
    %set(axs{iax},'xtickmode','manual','xtick', []);
   
end

fig_m = figure('Name', ['MCMC E-Step: trend of means']);
for iax = 1:N_CUST
    axs{iax} = subplot(7, 1, iax);
    
    set(axs{iax}, 'fontsize', 7);
    set(axs{iax},'XTick',[]);
    set(axs{iax},'ylim',[0 340], 'ytick', 0:100:300);
    part_dat = dat(:, iax:N_CUST:end);
    for iline = 1:N_STEP+1
        n_sample = pop_size/N_STEP*iline;
        if iline == N_STEP+1
            %n_sample = 1;
            continue;
        end
        line(1:N_TIME, mean(part_dat(1:n_sample, :), 1)', ...
            'color', color{iline});
    end
end
legend('show');


tx = 0:24:168;
set(axs{N_CUST}, 'xtick', tx, 'xticklabel', cellstr(num2str(tx')));
    
    


end

