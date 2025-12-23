xyz_ind_gslib=discretize(input_gslib,Ncat);
if loglag
    lagsteplog=logspace(0,log10((max([NJ,NI,NK])-1)/2),lagsteplogspace);
    Lagsin=unique(floor(lagsteplog));
else
    Lagsin=1:lagstep:(max([NJ,NI,NK])-1);
end
for nd =1:NUMDIR
    textdisp= ['strike = ',num2str(alpha(nd)),';  dip = ',num2str(beta(nd))];
    disp(textdisp)
    [H0,HL_d]=entrogram_v2(xyz_ind_gslib,Lagsin,bandwidth_X,bandwidth_Y,bandwidth_Z,...
        alpha(nd),beta(nd),NI,NJ,NK,Ncat,NMC);
    eval(strcat('HL_d',num2str(nd),' = HL_d;'));
    eval(strcat('H0_d',num2str(nd),' = H0;'));
    dist_d = linspace(1,sqrt((bandwidth_X)^2+(bandwidth_Y)^2+(bandwidth_Y)^2),numel(HL_d));

    eval(strcat('dist_d',num2str(nd),' = dist_d;'));
end

col = ['r';'b';'g']; %colours for plotting entrograms
mark= ['o';'s';'x'];
figure('color','w')
subplot(2,1,1) % Plot the H0-normalised entrogram
hold on
for nd =1:NUMDIR
    eval(['HL = ',strcat('HL_d',num2str(nd)),';']);
    eval(['H0 = ',strcat('H0_d',num2str(nd)),';']);
    eval(['dist = ',strcat('dist_d',num2str(nd)),';']);
    plot1=plot([0;dist'/dx],[0;HL/H0],'Marker',mark(nd),'linewidth',1.,'markersize',5,'color',col(nd));
    xlabel('Lags');
    ylabel('H_R=H_L/H_0');
    h = legend('show');
    set(h,'location','southeast','fontsize',10, 'fontname', 'Arial')
    set(gca, 'fontsize',10,'linewidth',1,'FontName','Arial','box', 'off',...
        'tickdir','out','ticklength',[0.009,0.009]);
    H0_nd(nd)=H0;

    [A,B]=sort(abs(dist-numel(HL)));

    Hscale(nd) = trapz([0;dist(1:B(1))'./dx],[1;1-HL(1:B(1))/H0]);       % Calculate entropic scale
    normHscale(nd)=Hscale(nd)/(dist(end)./dx);            % Calculate normalized entropic scale
    if nHs_flag
        name2plot=horzcat('dir ',num2str(nd),' - nH_{S}= ',num2str(normHscale(nd),'%3.2f'));
    else
        name2plot=horzcat('dir ',num2str(nd),' - H_{S}= ',num2str(Hscale(nd),'%3.2f'));
    end
    set(plot1,'DisplayName',name2plot);
    ylim([0.0 1.1])
end

%set the dataset name and report relative and global entropy
tit=title(strcat('H_{R0} =',num2str(HL(1)/H0,'%3.2f'),' - H_0=',num2str(H0,'%3.2f')));
set(tit,'FontSize',12);

subplot(2,1,2) % Plot the not-normalised entrogram
hold on
for nd =1:NUMDIR
    eval([ 'HL= ',strcat('HL_d',num2str(nd)),';']);
    eval([ 'dist= ',strcat('dist_d',num2str(nd)),';']);

    plot1=plot([0;dist'/dx],[0;HL],'Marker',mark(nd),'linewidth',1.,'markersize',5,'color',col(nd));
    xlabel('Lags (Distance)');
    ylabel('H_L');

    name2plot=horzcat('H_0(dir ',num2str(nd),')\rightarrow',num2str(mean(HL(end-3:end)),'%3.2f'));
    set(plot1,'DisplayName',name2plot);
    h = legend('show');
    set(h,'location','southeast','fontsize',10, 'fontname', 'Arial')
    set(gca, "fontsize",10,'linewidth',1,'FontName','Arial','box', 'off',...
        'tickdir','out','ticklength',[0.009,0.009]);
end

str1=a(j).name;



