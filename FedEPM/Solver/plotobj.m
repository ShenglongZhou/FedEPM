function  plotobj(obj)
    figure('Renderer', 'painters', 'Position',[1100 400 370 320]);
    axes('Position', [0.07 0.15 0.83 0.8] ); 
    h1 = plot(1:length(obj),obj(1:end)); hold on
    %h2 = plot(1:length(obj2),obj2(1:end)); hold on
    grid on
    h1.LineWidth  = 1.5;   %h2.LineWidth  = 1.5;       
    h1.LineStyle  = '-';   %h2.LineStyle  = '--'; 
    h1.Color = '#3caea3';  %h2.Color = '#ed553b';
    axis([0 length(obj) 0.54 0.71]); grid on
    set(gca, 'YTick',0.55:0.05:0.7,'FontSize',12) ;
    yticklabels({'0.55','0.60','0.65','0.70'});   
    ly = ylabel('f(w)/m'); ly.Position(1) = ly.Position(1)+ 0.55;
    ax = gca; ax.YAxisLocation = 'right';
    xlabel('CR'); 
end
