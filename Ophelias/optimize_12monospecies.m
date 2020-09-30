function MSEV = optimize_12monospecies(Po,Gx)
    
    %orgorder = {'BH','CA','BU','PC','BO','BV','BT','EL','FP','CH','DP','ER'};
    
    G1 = Gx(1,:);
    G2 = Gx(2,:);
    G3 = Gx(3,:);
    G4 = Gx(4,:);
    G5 = Gx(5,:);
    G6 = Gx(6,:);
    G7 = Gx(7,:);
    G8 = Gx(8,:);
    G9 = Gx(9,:);
    G10 = Gx(10,:);
    G11 = Gx(11,:);
    G12 = Gx(12,:);
    
    timem1 = ((0:1:size(G1,2)-1)*30)./60;
    timem2 = ((0:1:size(G2,2)-1)*30)./60;
    timem3 = ((0:1:size(G3,2)-1)*30)./60;
    timem4 = ((0:1:size(G4,2)-1)*30)./60;
    timem5 = ((0:1:size(G5,2)-1)*30)./60;
    timem6 = ((0:1:size(G6,2)-1)*30)./60;
    timem7 = ((0:1:size(G7,2)-1)*30)./60;
    timem8 = ((0:1:size(G8,2)-1)*30)./60;
    timem9 = ((0:1:size(G9,2)-1)*30)./60;
    timem10 = ((0:1:size(G10,2)-1)*30)./60;
    timem11 = ((0:1:size(G11,2)-1)*30)./60;
    timem12 = ((0:1:size(G12,2)-1)*30)./60;
    
    icx = zeros(1,12);
    %MONOCULTURE STRAIN 1
    if G1(1)==0
       xnot = G1(2);
    else 
       xnot = G1(1);
    end
    ic1 = icx;
    ic1(1) = xnot;
    MoutS1=LV_12member(Po(:),ic1,max(timem1));
    for i = 1:length(timem1)
        S1(i) = interp1(MoutS1(:,1),MoutS1(:,2),timem1(i));
    end   
    MSEVS1 = immse(S1,G1);
    
%     figure
%     plot(timem1,S1,'ro-','LineWidth',2);
%     hold on
%     plot(timem1,G1,'bo-','LineWidth',2);
%     set(gca,'FontSize',20);
%     title('G1');
         
    %MONOCULTURE STRAIN 2
    if G2(1)==0
       xnot = G2(2);
    else
       xnot = G2(1);
    end
    ic2 = icx;
    ic2(2) = xnot;
    MoutS2=LV_12member(Po(:),ic2,max(timem2));
    for i = 1:length(timem2)
        S2(i) = interp1(MoutS2(:,1),MoutS2(:,3),timem2(i));
    end
    MSEVS2 = immse(S2,G2);
    
%     figure
%     plot(timem1,S2,'ro-','LineWidth',2);
%     hold on
%     plot(timem1,G2,'bo-','LineWidth',2);
%     set(gca,'FontSize',20);
%     title('G2');
    
    %MONOCULTURE STRAIN 3
    if G3(1)==0
       xnot = G3(2);
    else
       xnot = G3(1);
    end
    ic3 = icx;
    ic3(3) = xnot;
    MoutS3=LV_12member(Po(:),ic3,max(timem3));
    for i = 1:length(timem3)
        S3(i) = interp1(MoutS3(:,1),MoutS3(:,4),timem3(i));
    end
    MSEVS3 = immse(S3,G3);
    
%     figure
%     plot(timem1,S3,'ro-','LineWidth',2);
%     hold on
%     plot(timem1,G3,'bo-','LineWidth',2);
%     set(gca,'FontSize',20);
%     title('G3');
    
    %MONOCULTURE STRAIN 4
    if G4(1)==0
       xnot = G4(2);
    else
       xnot = G4(1);
    end
    ic4 = icx;
    ic4(4) = xnot;
    MoutS4=LV_12member(Po(:),ic4,max(timem4));
    for i = 1:length(timem4)
        S4(i) = interp1(MoutS4(:,1),MoutS4(:,5),timem4(i));
    end
    MSEVS4 = immse(S4,G4);
    
%     figure
%     plot(timem1,S4,'ro-','LineWidth',2);
%     hold on
%     plot(timem1,G4,'bo-','LineWidth',2);
%     set(gca,'FontSize',20);
%     title('G4');
    
    %MONOCULTURE STRAIN 5
    if G5(1)==0
       xnot = G5(2);
    else
       xnot = G5(1);
    end
    ic5 = icx;
    ic5(5) = xnot;
    MoutS5=LV_12member(Po(:),ic5,max(timem5));
    for i = 1:length(timem5)
        S5(i) = interp1(MoutS5(:,1),MoutS5(:,6),timem5(i));
    end
    MSEVS5 = immse(S5,G5);
    
%     figure
%     plot(timem1,S5,'ro-','LineWidth',2);
%     hold on
%     plot(timem1,G5,'bo-','LineWidth',2);
%     set(gca,'FontSize',20);
%     title('G5');
    
    
    %MONOCULTURE STRAIN 6
    if G6(1)==0
       xnot = G6(2);
    else
       xnot = G6(1);
    end
    ic6 = icx;
    ic6(6) = xnot;
    MoutS6=LV_12member(Po(:),ic6,max(timem6));
    for i = 1:length(timem6)
        S6(i) = interp1(MoutS6(:,1),MoutS6(:,7),timem6(i));
    end
    MSEVS6 = immse(S6,G6);
    
    %MONOCULTURE STRAIN 7
    if G7(1)==0
       xnot = G7(2);
    else
       xnot = G7(1);
    end
    ic7 = icx;
    ic7(7) = xnot;
    MoutS7=LV_12member(Po(:),ic7,max(timem7));
    for i = 1:length(timem7)
        S7(i) = interp1(MoutS7(:,1),MoutS7(:,8),timem7(i));
    end
    MSEVS7 = immse(S7,G7);
        
    %MONOCULTURE STRAIN 8
    if G8(1)==0
       xnot = G8(2);
    else
       xnot = G8(1);
    end
    ic8 = icx;
    ic8(8) = xnot;
    MoutS8=LV_12member(Po(:),ic8,max(timem8));
    for i = 1:length(timem8)
        S8(i) = interp1(MoutS8(:,1),MoutS8(:,9),timem8(i));
    end
    MSEVS8 = immse(S8,G8);
    
    %MONOCULTURE STRAIN 9
    if G9(1)==0
       xnot = G9(2);
    else
       xnot = G9(1);
    end
    ic9 = icx;
    ic9(9) = xnot;
    MoutS9=LV_12member(Po(:),ic9,max(timem9));
    for i = 1:length(timem9)
        S9(i) = interp1(MoutS9(:,1),MoutS9(:,10),timem9(i));
    end
    MSEVS9 = immse(S9,G9);
    
    %MONOCULTURE STRAIN 10
    if G10(1)==0
       xnot = G10(2);
    else
       xnot = G10(1);
    end
    ic10 = icx;
    ic10(10) = xnot;
    MoutS10=LV_12member(Po(:),ic10,max(timem10));
    for i = 1:length(timem10)
        S10(i) = interp1(MoutS10(:,1),MoutS10(:,11),timem10(i));
    end
    MSEVS10 = immse(S10,G10);
    
    %MONOCULTURE STRAIN 11
    if G11(1)==0
       xnot = G11(2);
    else
       xnot = G11(1);
    end
    ic11 = icx;
    ic11(11) = xnot;
    MoutS11=LV_12member(Po(:),ic11,max(timem11));
    for i = 1:length(timem11)
        S11(i) = interp1(MoutS11(:,1),MoutS11(:,12),timem11(i));
    end
    MSEVS11 = immse(S11,G11);
    
    %MONOCULTURE STRAIN 12
    if G12(1)==0
       xnot = G12(2);
    else
       xnot = G12(1);
    end
    ic12 = icx;
    ic12(12) = xnot;
    MoutS12=LV_12member(Po(:),ic12,max(timem12));
    for i = 1:length(timem12)
        S12(i) = interp1(MoutS12(:,1),MoutS12(:,13),timem12(i));
    end
    MSEVS12 = immse(S12,G12);
    
    orgorder = {'BH','CA','BU','PC','BO','BV','BT','EL','FP','CH','DP','ER'};
           
%     figure
%     plot(timem1,S6,'ro-','LineWidth',2);
%     hold on
%     plot(timem1,G6,'bo-','LineWidth',2);
%     set(gca,'FontSize',20);
%     title('G6');

    tx{1} = timem1;
    tx{2} = timem2;
    tx{3} = timem3;
    tx{4} = timem4;
    tx{5} = timem5;
    tx{6} = timem6;
    tx{7} = timem7;
    tx{8} = timem8;
    tx{9} = timem9;
    tx{10} = timem10;
    tx{11} = timem11;
    tx{12} = timem12;
    
    datav{1} = [S1; G1];
    datav{2} = [S2; G2];
    datav{3} = [S3; G3];
    datav{4} = [S4; G4];
    datav{5} = [S5; G5];
    datav{6} = [S6; G6];
    datav{7} = [S7; G7];
    datav{8} = [S8; G8];
    datav{9} = [S9; G9];
    datav{10} = [S10; G10];
    datav{11} = [S11; G11];
    datav{12} = [S12; G12];
       
%     figure;
%     for i = 1:12
%         subplot(3,4,i);
%         plot(tx{i},datav{i},'LineWidth',3)
%         hold on;
%         set(gca,'FontSize',14);
%         legend(orgorder{i});
%     end
    
    %MSEV = (MSEVS1+MSEVS2+MSEVS3+MSEVS4+MSEVS5+MSEVS6+MSEVS7+MSEVS8+MSEVS9+MSEVS10+MSEVS11+MSEVS12)./12;
    MSEV = MSEVS1+MSEVS2+MSEVS3+MSEVS4+MSEVS5+MSEVS6+MSEVS7+MSEVS8+MSEVS9+MSEVS10+MSEVS11+MSEVS12;
    
end
