% This script visualizes the results of statistical tests (V-test (mu=0)
% and Rayleigh test) on the behavior simulation (same test as the
% experimental paper) from server.
%
% Tianshu Li
% Dec. 6, 2022

fntsz = 18; % font size
stim_names = {'control0','DoG29','bar40','DoG69'};
legendnames = {'control', '29^\circ DoG', '40^\circ bar', '69^\circ DoG'}; % used for legend label
nstim=length(stim_names);
stats = struct('stim',stim_names); % save results in separate structures  
figure(1); clf; % box plot
figure(2); clf; % mean
for kname = 1:nstim
    % load data
    stim_name = stim_names{kname};
    load(sprintf('data/pvalues_Urchin_behavior_%s.mat',stim_name));

    % reshape the cell arrays to matrices
    p_vtest_mat = zeros(length(Nanimals),nsessions);
    p_Rayleigh_mat = p_vtest_mat; % 17x100: 17 is the number of sample sizes (40:10:200) and 100 is the number of urchins (=initial orientations) in each sample. 
    vtest_stat_mat = p_vtest_mat;
    rtest_stat_mat = p_vtest_mat;
    for ksession = 1:nsessions
        p_vtest_mat(:,ksession) = p_vtest{ksession}';
        p_Rayleigh_mat(:,ksession) = p_Rayleigh{ksession}';
        vtest_stat_mat(:,ksession) = vtest_stat{ksession}';
        rtest_stat_mat(:,ksession) = rtest_stat{ksession}';
    end

    % ===  generate plots  ===
    % ---  1. bar plot  ---
    figure(1);
    % V-test
    subplot(2,2,1); hold on;
    boxchart(p_vtest_mat');
    xticklabels(string(Nanimals));
    figset(gca,'n (sample size)','p-value',fntsz);
    title('V-test','FontSize',fntsz+2);

    subplot(2,2,3); hold on;
    boxchart(vtest_stat_mat');
    xticklabels(string(Nanimals));
    figset(gca,'n (sample size)','V-test statistic',fntsz);

    % Rayleigh test
    subplot(2,2,2); hold on;
    boxchart(p_Rayleigh_mat');
    xticklabels(string(Nanimals));
    figset(gca,'n (sample size)','p-value',fntsz);
    title('Rayleigh test','FontSize',fntsz+2);

    subplot(2,2,4); hold on;
    boxchart(rtest_stat_mat');
    xticklabels(string(Nanimals));
    figset(gca,'n (sample size)','Rayleigh test statistic',fntsz);

    % ---  2. mean  ---
    % average observables:
    obs1pval = mean(p_vtest_mat,2); % P-value of V test
    obs1sval = mean(vtest_stat_mat,2); % statistic of V test
    obs2pval = mean(p_Rayleigh_mat,2); % P value of Ra. test
    obs2sval = mean(rtest_stat_mat,2); % statistics of Ra. test
    
    figure(2);
    % V-test
    subplot(1,2,1); hold on;
    plot(Nanimals,obs1pval,'LineWidth',2);
    figset(gca,'n (sample size)','averaged p-value',fntsz);
    title('V-test','FontSize',fntsz+2);

    % Rayleigh test
    subplot(1,2,2); hold on;
    plot(Nanimals,obs2pval,'LineWidth',2);
    figset(gca,'n (sample size)','averaged p-value',fntsz);
    title('Rayleigh test','FontSize',fntsz+2);
    
    % save desired results for each stimulus 
    stats(kname).Vtest=[Nanimals; obs1sval' ; obs1pval']'; 
    stats(kname).Ratest=[Nanimals; obs2sval' ; obs2pval']';
end
figure(1);
legend(legendnames);
legend 'boxoff';

figure(2);
legend(legendnames);
legend 'boxoff';
for k = 1:2
    subplot(1,2,k);
    xlim([min(Nanimals),max(Nanimals)]);
    xticks(Nanimals);
    xticklabels(string(Nanimals));
end

%% print results to consolle
% print p-values to consolle for the last stimulus plotted:
% [Nanimals ; mean(p_vtest_mat,2)']'

% results for DoG 69deg stimulus (for n<=100)
idx=find(Nanimals==100); % find index of sample with 100 animals
idxname=find(strcmpi(stim_names,'DoG69')); % find index of DoG69 stimulus
disp('V test mean statistic and P values:');
stats(idxname).Vtest(1:idx,:)
disp('Rayleigh test mean statistic and P values:');
stats(idxname).Ratest(1:idx,:)

% p-values for all stimuli for n=100
disp('summary of P-values for all stimuli with n=100:');
for k=1:numel(stats)
    fprintf('%s: (V-test) %.3f, (Rayleigh test) %.3f\n',stats(k).stim,stats(k).Vtest(idx,end),stats(k).Ratest(idx,end));
end
















