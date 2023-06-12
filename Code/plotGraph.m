function [figure] = plotGraph(G,X,X_star, legend1, legend2, title_, save_, name)
n = size(X,1);
X_star = X_star - repmat(mean(X_star),n,1);
XX = X - repmat(mean(X),n,1);

if size(X,2) == 3

    figure  = plot(G,'XData',X_star(:,1),'YData',X_star(:,2),'ZData',X_star(:,3),'LineStyle','--','DisplayName',legend1)
    hold on
    plot(G,'XData',XX(:,1),'YData',XX(:,2),'ZData',XX(:,3),'LineStyle','-.','DisplayName',legend2)
    title(title)
    legend()
    hold off;

    if save_
        saveas(fig,name)
    end

elseif size(X,2) == 2
    
    figure = plot(G,'XData',X_star(:,1),'YData',X_star(:,2),'LineStyle','--','DisplayName',legend1)
    hold on
    plot(G,'XData',XX(:,1),'YData',XX(:,2),'LineStyle','-.','DisplayName',legend2)
    title(title_);
    legend();
    hold off;

    if save_
        saveas(figure,name)
    end
end
end