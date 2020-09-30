function cost = objectiveFunction_2(Po,lambda)
        lam = lambda;   
        
        orgorder = {'BH','CA','BU','PC','BO','BV','BT','EL','FP','CH','DP','ER'};
        %MONOCULTURE GROWTH DATA
        load growth_data_0.mat
        sg0 = singlegrowthn;

        for i = 1:length(orgorder)
            im(i) = find(ismember(singlen,orgorder{i}));
        end
        Gx = sg0(im,:);
                
        load PW_rel_abund_matrices_PW19_map.mat

        %MONOCULTURE 
        MSEVm = optimize_12monospecies(Po,Gx);
        
        %loop over all community measurements
        indx = 1:66; %JUST PAIRWISE
        for i = 1:length(indx)
            m = indx(i);
            
            if m <= 66 %Pairwise
                namex = uniquename{m};
                npair = m;
                nL1O6 = 0;
                nTM1 = 0;
                nL1O7 = 0;
                nPW22 = 0;
                if ~isnan(ABM(m,end))
                    MSEVc = optimize_12community(Po,TM(m,:),ABM(m,:),npair,namex);
                else
                    MSEVc = optimize_12community(Po,TM(m,1:end-1),ABM(m,1:end-1),npair,namex);
                end
            end
            
            MSEVt(i) = MSEVc;
        end

        out = sum(MSEVt) + MSEVm; %compute sum, all data sets equal
 
        Po = Po(:)';
        %%%REGULARIZATION
        %cost = out + lam*(sumabs(Po(:))); %L^1 regularization
        cost = out; %L^1 regularization
        disp(cost);
end
