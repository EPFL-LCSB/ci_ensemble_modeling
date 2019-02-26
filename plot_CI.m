function [] = plot_CI(means, CIs, plot_title)

[~,m_vec_order] = sort(abs(means),'descend');

Low_tail = CIs(1,:);
Upp_tail = CIs(2,:);

isSignif=[Upp_tail(m_vec_order)<0 | Low_tail(m_vec_order)>0];
figure
errorbar(find(~isSignif),means(m_vec_order(~isSignif)),means(m_vec_order(~isSignif))-Low_tail(m_vec_order(~isSignif))...
    ,Upp_tail(m_vec_order(~isSignif))-means(m_vec_order(~isSignif)),'db','MarkerFaceColor','b','MarkerSize',8,'LineWidth',1.5); 
hold on
errorbar(find(isSignif),means(m_vec_order(isSignif)),means(m_vec_order(isSignif))-Low_tail(m_vec_order(isSignif))...
    ,Upp_tail(m_vec_order(isSignif))-means(m_vec_order(isSignif)),'dr','MarkerFaceColor','r','MarkerSize',8,'LineWidth',1.5); 
xlabel('Variables')
ylabel('Means and CIs')
if nargin == 3
    title(plot_title)
end
set(gca,'fontsize',22)
set(gca,'fontweight','bold')
grid on

end