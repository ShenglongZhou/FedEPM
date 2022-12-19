function  plot4k0(out,k0,m)
    figure('Renderer', 'painters', 'Position',[1100 400 370 350]);
    axes('Position', [0.07 0.14 0.83 0.8] ); 
    styles  = {'-','-.',':','-.','--'}; 
    for i  = 1 : nnz(k0)
        y  = out{i}.objx;   
        pl = plot(log2(1:length(y)), y/m, styles{i}); hold on 
        pl.LineWidth = 2;
        leg{i} = strcat('k_0=',num2str(k0(i)));
    end 
    axis([0 5 0.54 0.71]); grid on
    set(gca, 'YTick',0.55:0.05:0.7,'FontSize',10) ;
    yticklabels({'0.55','0.60','0.65','0.70'});   
    ly = ylabel('f(w)/m'); ly.Position(1) = ly.Position(1)+0.2;
    set(gca, 'XTick',0:1:5,'FontSize',10) ;
    ax = gca; ax.YAxisLocation = 'right';
    xticklabels({'2^0','2^1','2^2','2^3','2^4','2^5'}); xlabel('CR');   
    title(strcat('m = ',num2str(m)),'FontWeight','normal') 
    legend(leg,'Location','NorthEast');
 
end
