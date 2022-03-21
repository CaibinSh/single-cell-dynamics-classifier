function plot2p(p53,p21,Xaxis,varargin)
%% PLOT2P function is designed to plot 25th percentile, 75th percentile and median level of both p53 and p21 in a single plot
% written by Caibin Sheng, Loewer lab, TU Darmstadt

    if nargin==3
        druginsertion=1; % default druginsertion is the first time point
        
    elseif nargin==4
        druginsertion=varargin{1};
    end
        druginsertiontime = Xaxis(druginsertion);
        
        temp2_p53 = prctile(p53,[25 50 75],2);
        temp2_p21 = prctile(p21,[25 50 75],2);
           
        [AX,hLine1,hLine2] = plotyy(Xaxis,temp2_p53(:,2),Xaxis,temp2_p21(:,2));set(hLine1,'Color',[0 0.57 0.57]);set(hLine1,'LineWidth',1.8);set(hLine2,'Color',[1 0.43 0.71]);set(hLine2,'LineWidth',1.8);
        set(AX(2),'YLim',[0 1000]),set(AX(2),'XLim',[0 max(Xaxis)]);
        set(AX,{'ycolor'},{[0 0.57 0.57];[1 0.43 0.71]})
        
        hold(AX(2)), plot(AX(2),Xaxis,temp2_p21(:,1),'color',[1 0.43 0.71])
        hold on, plot(AX(2),Xaxis,temp2_p21(:,3),'color',[1 0.43 0.71])
        set(AX(2),'YTick',[0 500 1000])        
        ylabel(AX(2),'p21 level (a.u.)')
               
        hold off
        
        hold(AX(1)),plot(AX(1),Xaxis,temp2_p53(:,1),'color',[0 0.57 0.57])
        hold on, plot(AX(1),Xaxis,temp2_p53(:,3),'color',[0 0.57 0.57])
        set(AX(1),'YLim',[0 1000]) 
        
        linkaxes([AX(2),AX(1)],'x')
        set(AX(1),'YTick',[0 500 1000])
        text(3,800,['cell #: ', num2str(size(p21,2))],'FontSize',8)
        
        ylabel(AX(1),'p53 level (a.u.)')
          
        hold on
        if druginsertiontime>1
            line([druginsertiontime druginsertiontime],[0 1000])
        end
        ax = AX(1);
        a = sort(unique([fliplr(Xaxis(druginsertion):-10:0),Xaxis(druginsertion):10:max(Xaxis)]));
        ax.XTick = a;
        ax.XTickLabel = a-druginsertiontime;
        box on
        xlabel('time after irradiation (h)')
end