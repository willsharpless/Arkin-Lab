function MSEV = optimize_12community(Po,timeseq,relabund,npair,commname)

        orgorder = {'BH','CA','BU','PC','BO','BV','BT','EL','FP','CH','DP','ER'};
        
        %LOAD PW DATA
        load PW_rel_abund_matrices_PW19_map.mat
        load growth_data_0.mat
        pg0 = pairgrowthn;
        load growth_data_1.mat
        pg2 = pairgrowthn;
        load growth_data_2.mat
        pg4 = pairgrowthn;
        
        ict = zeros(1,12);
        t1 = timeseq(1:3);
        t2 = timeseq(3:5)-timeseq(3);
        t3 = timeseq(5:end)-timeseq(5);
        
        %PAIRWISE
        if npair > 0
           a0 = interpolate_monoculture_data(pg0(npair,:),t1(1:3),(0:1:size(pg0,2)-1)*30./60);
           a2 = interpolate_monoculture_data(pg2(npair,:),t2(2:end),(0:1:size(pg2,2)-1)*30./60);
           a4 = interpolate_monoculture_data(pg4(npair,:),t3(2:end),(0:1:size(pg4,2)-1)*30./60);

            o1 = find(ismember(orgorder,commname(1:2)));
            o2 = find(ismember(orgorder,commname(3:4)));

            ict(o1) = relabund(1);
            ict(o2) = 1-relabund(1);
            ncells = 0.02;

            relabundmat = zeros(length(relabund),12);
            for i = 1:length(relabund)
                relabundmat(i,o1) = relabund(i);
                relabundmat(i,o2) = 1-relabund(i);
            end

        end

        check = isnan(a0);
        ix = find(check==1);
        if ~isempty(ix)
           a0(ix) = 0;
        end

        NG0 = a0;
        NG2 = [a0(end)./20 a2];
        NG4 = [a2(end)./20 a4];

        nn = [NG0(1:3) NG2(2:end) NG4(2:end)];

        for i = 1:size(nn,2)
            absabund(i,:) = relabundmat(i,:).*nn(i);
        end


        %TRANSFER 0
        ic0 = ncells*ict;
        Tmax = t1(end);
        Mout0=LV_12member(Po(:),ic0,Tmax);
        ta = Mout0(:,1);
        x1oa = Mout0(:,2);
        x2oa = Mout0(:,3);
        x3oa = Mout0(:,4);
        x4oa = Mout0(:,5);
        x5oa = Mout0(:,6);
        x6oa = Mout0(:,7);
        x7oa = Mout0(:,8);
        x8oa = Mout0(:,9);
        x9oa = Mout0(:,10);
        x10oa = Mout0(:,11);
        x11oa = Mout0(:,12);
        x12oa = Mout0(:,13);


        for i = 1:length(t1)
            ao1(i) = interp1(ta,x1oa,t1(i));
            ao2(i) = interp1(ta,x2oa,t1(i));
            ao3(i) = interp1(ta,x3oa,t1(i));
            ao4(i) = interp1(ta,x4oa,t1(i));
            ao5(i) = interp1(ta,x5oa,t1(i));
            ao6(i) = interp1(ta,x6oa,t1(i));
            ao7(i) = interp1(ta,x7oa,t1(i));
            ao8(i) = interp1(ta,x8oa,t1(i));
            ao9(i) = interp1(ta,x9oa,t1(i));
            ao10(i) = interp1(ta,x10oa,t1(i));
            ao11(i) = interp1(ta,x11oa,t1(i));
            ao12(i) = interp1(ta,x12oa,t1(i));
        end

        %TRANSFER 2
        ic2 = [ao1(end) ao2(end) ao3(end) ao4(end) ao5(end) ao6(end) ao7(end) ao8(end) ao9(end) ao10(end) ao11(end) ao12(end)]./20;
        Tmax = t2(end);
        Mout1=LV_12member(Po(:),ic2,Tmax);
        tb = Mout1(:,1);
        x1ob=Mout1(:,2);
        x2ob=Mout1(:,3);
        x3ob=Mout1(:,4);
        x4ob=Mout1(:,5);
        x5ob=Mout1(:,6);
        x6ob=Mout1(:,7);
        x7ob=Mout1(:,8);
        x8ob=Mout1(:,9);
        x9ob=Mout1(:,10);
        x10ob=Mout1(:,11);
        x11ob=Mout1(:,12);
        x12ob=Mout1(:,13);

        for i = 1:length(t2)
            a21(i) = interp1(tb,x1ob,t2(i));
            a22(i) = interp1(tb,x2ob,t2(i));
            a23(i) = interp1(tb,x3ob,t2(i));
            a24(i) = interp1(tb,x4ob,t2(i));
            a25(i) = interp1(tb,x5ob,t2(i));
            a26(i) = interp1(tb,x6ob,t2(i));
            a27(i) = interp1(tb,x7ob,t2(i));
            a28(i) = interp1(tb,x8ob,t2(i));
            a29(i) = interp1(tb,x9ob,t2(i));
            a210(i) = interp1(tb,x10ob,t2(i));
            a211(i) = interp1(tb,x11ob,t2(i));
            a212(i) = interp1(tb,x12ob,t2(i));
        end


        %TRANSFER 4
        ic4 = [a21(end) a22(end) a23(end) a24(end) a25(end) a26(end) a27(end) a28(end) a29(end) a210(end) a211(end) a212(end)]./20;
        Tmax = t3(end);
        Mout2=LV_12member(Po(:),ic4,Tmax);
        tc = Mout2(:,1);
        x1oc=Mout2(:,2);
        x2oc=Mout2(:,3);
        x3oc=Mout2(:,4);
        x4oc=Mout2(:,5);
        x5oc=Mout2(:,6);
        x6oc=Mout2(:,7);
        x7oc=Mout2(:,8);
        x8oc=Mout2(:,9);
        x9oc=Mout2(:,10);
        x10oc=Mout2(:,11);
        x11oc=Mout2(:,12);
        x12oc=Mout2(:,13);

        for i = 1:length(t3)
            a41(i) = interp1(tc,x1oc,t3(i));
            a42(i) = interp1(tc,x2oc,t3(i));
            a43(i) = interp1(tc,x3oc,t3(i));
            a44(i) = interp1(tc,x4oc,t3(i));
            a45(i) = interp1(tc,x5oc,t3(i));
            a46(i) = interp1(tc,x6oc,t3(i));
            a47(i) = interp1(tc,x7oc,t3(i));
            a48(i) = interp1(tc,x8oc,t3(i));
            a49(i) = interp1(tc,x9oc,t3(i));
            a410(i) = interp1(tc,x10oc,t3(i));
            a411(i) = interp1(tc,x11oc,t3(i));
            a412(i) = interp1(tc,x12oc,t3(i));
        end


        X1v = [ao1 a21(2:end) a41(2:end)];
        X2v = [ao2 a22(2:end) a42(2:end)];
        X3v = [ao3 a23(2:end) a43(2:end)];
        X4v = [ao4 a24(2:end) a44(2:end)];
        X5v = [ao5 a25(2:end) a45(2:end)];
        X6v = [ao6 a26(2:end) a46(2:end)];
        X7v = [ao7 a27(2:end) a47(2:end)];
        X8v = [ao8 a28(2:end) a48(2:end)];
        X9v = [ao9 a29(2:end) a49(2:end)];
        X10v = [ao10 a210(2:end) a410(2:end)];
        X11v = [ao11 a211(2:end) a411(2:end)];
        X12v = [ao12 a212(2:end) a412(2:end)];

        Idata = [X1v; X2v; X3v; X4v; X5v; X6v; X7v; X8v; X9v; X10v; X11v; X12v];


        for i = 1:size(Idata,2)
            Idatarel(:,i) = Idata(:,i)./sum(Idata(:,i));
        end


        %%ABSOLUTE ABUNDANCE PLOT
        %figure(1);
%         figure;
%         for i = 1:12
%             subplot(4,3,i)
%             plot(timeseq,Idata(i,:),'o-','LineWidth',2,'MarkerSize',10);
%             hold on
%             plot(timeseq,absabund(:,i),'s-','LineWidth',2,'MarkerSize',10);
%             hold on;
%             if i == 1
%                 legend('model','data');
%             end
%             set(gca,'FontSize',16);
%             title(orgorder{i});
%         end
        
        %%RELATIVE ABUNDANCE PLOT
%         figure;
%         for i = 1:12
%             subplot(4,3,i)
%             plot(timeseq,Idatarel(i,:),'bo-','LineWidth',2,'MarkerSize',10);
%             hold on
%             plot(timeseq,relabundmat(:,i),'rs-','LineWidth',2,'MarkerSize',10);
%             hold on;
%             if i == 1
%                 legend('model','data');
%             end
%             set(gca,'FontSize',16);
%             title(orgorder{i});
%         end

        MSEV1 = nanmean((absabund(:,1)-X1v').^2);
        MSEV2 = nanmean((absabund(:,2)-X2v').^2);
        MSEV3 = nanmean((absabund(:,3)-X3v').^2);
        MSEV4 = nanmean((absabund(:,4)-X4v').^2);
        MSEV5 = nanmean((absabund(:,5)-X5v').^2);
        MSEV6 = nanmean((absabund(:,6)-X6v').^2);
        MSEV7 = nanmean((absabund(:,7)-X7v').^2);
        MSEV8 = nanmean((absabund(:,8)-X8v').^2);
        MSEV9 = nanmean((absabund(:,9)-X9v').^2);
        MSEV10 = nanmean((absabund(:,10)-X10v').^2);
        MSEV11 = nanmean((absabund(:,11)-X11v').^2);
        MSEV12 = nanmean((absabund(:,12)-X12v').^2);

        MSEV = MSEV1 + MSEV2 + MSEV3 + MSEV4 + MSEV5 + MSEV6 + MSEV7 + MSEV8 + MSEV9 + MSEV10 + MSEV11 + MSEV12;
end
