 function  plotobj(obj)
    figure('Renderer', 'painters', 'Position',[1100 400 370 320]);
    axes('Position', [0.07 0.15 0.83 0.8] ); 
    h1 = plot(1:length(obj),obj(1:end)); 
    hold on; grid on; xlabel('CR','FontSize',11)
    h1.LineWidth  = 1.5;     
    h1.LineStyle  = '-';  
    h1.Color      = '#3caea3';  
    axis([0 length(obj) 0.54 0.71]); grid on
    set(gca, 'YTick',0.55:0.05:0.7,'FontSize',11) ;
    yticklabels({'0.55','0.60','0.65','0.70'});   
    ly = ylabel('f(w)/m','FontSize',11); 
    ly.Position(1) = ly.Position(1)+ 0.45;
    ax = gca; ax.YAxisLocation = 'right';
    legend('FedEPM')
end
