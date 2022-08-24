function  plot4k0(out,k0,m)
    figure('Renderer', 'painters', 'Position',[900 200 497.5 460]);
    axes('Position', [0.06 0.11 0.86 0.85] ); 
    styles  = {'-','-.',':','-.','--'}; 
    for i  = 1 : nnz(k0)
        y  = out{i}.objx;   
        pl = plot(log2(1:length(y)), y/m, styles{i}); hold on 
        pl.LineWidth = 2;
    end 
    axis([0 5 0.54 0.71]); grid on
    set(gca, 'YTick',0.55:0.05:0.7,'FontSize',12) ;
    yticklabels({'0.55','0.60','0.65','0.70'});   
    ly = ylabel('f(w)/m'); ly.Position(1) = ly.Position(1)+0.1;
    set(gca, 'XTick',0:1:5,'FontSize',12) ;
    ax = gca; ax.YAxisLocation = 'right';
    xticklabels({'2^0','2^1','2^2','2^3','2^4','2^5'}); 
    xlabel('CR');   
    title(strcat('m = ',num2str(m)),'FontWeight','normal') 
    legend('k_0=4','k_0=8','k_0=12','k_0=16','k_0=20','Location','NorthEast');
 
end
