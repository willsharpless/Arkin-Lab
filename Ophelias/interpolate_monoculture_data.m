function vals = interpolate_monoculture_data(strain1,tp,tx)
    
    if max(tp)>max(tx)
       tp(end) = tx(end);
    end

    org1=strain1;
    for i = 1:length(tp)
    xx(i) = interp1(tx,org1,tp(i));
    end
    
    if isnan(xx(end))
       iltp = max(find(~isnan(strain1)));
       xx(end) = strain1(iltp);
    end 

    vals = xx;
    
%     figure;
%     plot(tx,org1,'Color',[rand rand rand],'LineWidth',3);
%     hold on;
%     for i = 1:length(vals)
%        plot([tp(i) tp(i)],[vals(i) vals(i)],'k*','MarkerSize',14);
%        hold on;
%     end

end
