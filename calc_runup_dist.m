function [R2]=calc_runup_dist(runup)
% This script calculates the runup distribution 
%runup_tmsrs, matrix containing all timeserires to be analyzed.

i=1;
%L = filter(ones(3,1)/3,1,[runs(k).one(:,i)' zeros(1,1)']);
%out1 = L(2:end);

[pks, locs]  = peakseek(runup, 4, 0);  % peak thresholds as deterimined by Pat
mx(i)=max(pks);
avg(i)=mean(pks);
mn(i)=min(pks);
sd(i)=std(pks);
%%Filter before getting peaks
% L = filter(ones(3,1)/3,1,[runup_tmsrs(:,i)' zeros(1,1)']);
% out1 = L(2:end);
% [pks,locs]  = findpeaks( out1);
% mx2=max(pks);

% %Plot Runup Exceedance Distribution
% figure(1)
% xlabel('Runup (m)','FontSize',20,'FontWeight','bold');
% ylabel('P(x)','FontSize',20,'FontWeight','bold');
% set(gca,'FontSize',20)
% %title(['Cumulative Distributions, H_s= ' num2str(Hs1(i)) 'm  T_p= ' num2str(Tp1(i)) ' sec'])
% box on
% grid on

asort=sort(pks);
nn=length(pks);
rank=1:length(pks);
for ii=1:length(pks)
    exprob(ii)=rank(ii)/nn;
end
% hold on
% plot(asort,exprob,'k--','LineWidth',0.5);

R2(i) = interp1(exprob, asort, 0.98);

    
 




