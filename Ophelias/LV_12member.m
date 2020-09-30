function Mout = LV_12member(params,icx,Tmax,orgs)

    r1 = params(1);
    a1_1 = params(2);
    a1_2 = params(3);
    a1_3 = params(4);
    a1_4 = params(5);
    a1_5 = params(6);
    a1_6 = params(7);
    a1_7 = params(8);
    a1_8 = params(9);
    a1_9 = params(10);
    a1_10 = params(11);
    a1_11 = params(12);
    a1_12 = params(13);

    r2 = params(14);
    a2_1 = params(15);
    a2_2 = params(16);
    a2_3 = params(17);
    a2_4 = params(18);
    a2_5 = params(19);
    a2_6 = params(20);
    a2_7 = params(21);
    a2_8 = params(22);
    a2_9 = params(23);
    a2_10 = params(24);
    a2_11 = params(25);
    a2_12 = params(26);

    r3 = params(27);
    a3_1 = params(28);
    a3_2 = params(29);
    a3_3 = params(30);
    a3_4 = params(31);
    a3_5 = params(32);
    a3_6 = params(33);
    a3_7 = params(34);
    a3_8 = params(35);
    a3_9 = params(36);
    a3_10 = params(37);
    a3_11 = params(38);
    a3_12 = params(39);

    r4 = params(40);
    a4_1 = params(41);
    a4_2 = params(42);
    a4_3 = params(43);
    a4_4 = params(44);
    a4_5 = params(45);
    a4_6 = params(46);
    a4_7 = params(47);
    a4_8 = params(48);
    a4_9 = params(49);
    a4_10 = params(50);
    a4_11 = params(51);
    a4_12 = params(52);

    r5 = params(53);
    a5_1 = params(54);
    a5_2 = params(55);
    a5_3 = params(56);
    a5_4 = params(57);
    a5_5 = params(58);
    a5_6 = params(59);
    a5_7 = params(60);
    a5_8 = params(61);
    a5_9 = params(62);
    a5_10 = params(63);
    a5_11 = params(64);
    a5_12 = params(65);

    r6 = params(66);
    a6_1 = params(67);
    a6_2 = params(68);
    a6_3 = params(69);
    a6_4 = params(70);
    a6_5 = params(71);
    a6_6 = params(72);
    a6_7 = params(73);
    a6_8 = params(74);
    a6_9 = params(75);
    a6_10 = params(76);
    a6_11 = params(77);
    a6_12 = params(78);

    r7 = params(79);
    a7_1 = params(80);
    a7_2 = params(81);
    a7_3 = params(82);
    a7_4 = params(83);
    a7_5 = params(84);
    a7_6 = params(85);
    a7_7 = params(86);
    a7_8 = params(87);
    a7_9 = params(88);
    a7_10 = params(89);
    a7_11 = params(90);
    a7_12 = params(91);

    r8 = params(92);
    a8_1 = params(93);
    a8_2 = params(94);
    a8_3 = params(95);
    a8_4 = params(96);
    a8_5 = params(97);
    a8_6 = params(98);
    a8_7 = params(99);
    a8_8 = params(100);
    a8_9 = params(101);
    a8_10 = params(102);
    a8_11 = params(103);
    a8_12 = params(104);

    r9 = params(105);
    a9_1 = params(106);
    a9_2 = params(107);
    a9_3 = params(108);
    a9_4 = params(109);
    a9_5 = params(110);
    a9_6 = params(111);
    a9_7 = params(112);
    a9_8 = params(113);
    a9_9 = params(114);
    a9_10 = params(115);
    a9_11 = params(116);
    a9_12 = params(117);

    r10 = params(118);
    a10_1 = params(119);
    a10_2 = params(120);
    a10_3 = params(121);
    a10_4 = params(122);
    a10_5 = params(123);
    a10_6 = params(124);
    a10_7 = params(125);
    a10_8 = params(126);
    a10_9 = params(127);
    a10_10 = params(128);
    a10_11 = params(129);
    a10_12 = params(130);

    r11 = params(131);
    a11_1 = params(132);
    a11_2 = params(133);
    a11_3 = params(134);
    a11_4 = params(135);
    a11_5 = params(136);
    a11_6 = params(137);
    a11_7 = params(138);
    a11_8 = params(139);
    a11_9 = params(140);
    a11_10 = params(141);
    a11_11 = params(142);
    a11_12 = params(143);

    r12 = params(144);
    a12_1 = params(145);
    a12_2 = params(146);
    a12_3 = params(147);
    a12_4 = params(148);
    a12_5 = params(149);
    a12_6 = params(150);
    a12_7 = params(151);
    a12_8 = params(152);
    a12_9 = params(153);
    a12_10 = params(154);
    a12_11 = params(155);
    a12_12 = params(156);

    ic1 = icx(1);
    ic2 = icx(2);
    ic3 = icx(3);
    ic4 = icx(4); 
    ic5 = icx(5);
    ic6 = icx(6);
    ic7 = icx(7);
    ic8 = icx(8);
    ic9 = icx(9);
    ic10 = icx(10);
    ic11 = icx(11);
    ic12 = icx(12);
    
    [t,y] = ode45(@dydt, [0 Tmax],[ic1 ic2 ic3 ic4 ic5 ic6 ic7 ic8 ic9 ic10 ic11 ic12]);
    
    for i = 1:size(y,1)
        x1n(i,:) = y(i,1)./(y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5)+y(i,6)+y(i,7)+y(i,8)+y(i,9)+y(i,10)+y(i,11)+y(i,12));
        x2n(i,:) = y(i,2)./(y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5)+y(i,6)+y(i,7)+y(i,8)+y(i,9)+y(i,10)+y(i,11)+y(i,12));
        x3n(i,:) = y(i,3)./(y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5)+y(i,6)+y(i,7)+y(i,8)+y(i,9)+y(i,10)+y(i,11)+y(i,12));
        x4n(i,:) = y(i,4)./(y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5)+y(i,6)+y(i,7)+y(i,8)+y(i,9)+y(i,10)+y(i,11)+y(i,12));
        x5n(i,:) = y(i,5)./(y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5)+y(i,6)+y(i,7)+y(i,8)+y(i,9)+y(i,10)+y(i,11)+y(i,12));
        x6n(i,:) = y(i,6)./(y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5)+y(i,6)+y(i,7)+y(i,8)+y(i,9)+y(i,10)+y(i,11)+y(i,12));
        x7n(i,:) = y(i,7)./(y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5)+y(i,6)+y(i,7)+y(i,8)+y(i,9)+y(i,10)+y(i,11)+y(i,12));
        x8n(i,:) = y(i,8)./(y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5)+y(i,6)+y(i,7)+y(i,8)+y(i,9)+y(i,10)+y(i,11)+y(i,12));
        x9n(i,:) = y(i,9)./(y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5)+y(i,6)+y(i,7)+y(i,8)+y(i,9)+y(i,10)+y(i,11)+y(i,12));
        x10n(i,:) = y(i,10)./(y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5)+y(i,6)+y(i,7)+y(i,8)+y(i,9)+y(i,10)+y(i,11)+y(i,12));
        x11n(i,:) = y(i,11)./(y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5)+y(i,6)+y(i,7)+y(i,8)+y(i,9)+y(i,10)+y(i,11)+y(i,12));
        x12n(i,:) = y(i,12)./(y(i,1)+y(i,2)+y(i,3)+y(i,4)+y(i,5)+y(i,6)+y(i,7)+y(i,8)+y(i,9)+y(i,10)+y(i,11)+y(i,12));   
    end
    
    %figure
    %plot(t,x1n,'LineWidth',3);
    %hold on
    %plot(t,x2n,'r-','LineWidth',3);
    %legend('x1','x2');
   
    %figure
    %plot(t,x1o,'LineWidth',3);
    %hold on
    %plot(t,x2o,'r-','LineWidth',3);
    %legend('x1','x2');
    
    Mout = [t y];
    
    %figure
    %plot(t,y,'LineWidth',4);
    %set(gca,'FontSize',20);
    %legend(orgs);
    %legend('x1','x2','x3','x4','x5','x6','x7','x8','x9','x10','x11','x12');
    
    function output = dydt(t,y)

        x1 = y(1);
        x2 = y(2);
        x3 = y(3);
        x4 = y(4);
        x5 = y(5);
        x6 = y(6);
        x7 = y(7);
        x8 = y(8);
        x9 = y(9);
        x10 = y(10);
        x11 = y(11);
        x12 = y(12);

        dydt(1) = r1*x1 + a1_1*x1^2 + a1_2*x1*x2 + a1_3*x1*x3 + a1_4*x1*x4 + a1_5*x1*x5 + a1_6*x1*x6 + a1_7*x1*x7 + a1_8*x1*x8 + a1_9*x1*x9 + a1_10*x1*x10 + a1_11*x1*x11 + a1_12*x1*x12;
        dydt(2) = r2*x2 + a2_2*x2^2 + a2_1*x1*x2 + a2_3*x2*x3 + a2_4*x2*x4 + a2_5*x2*x5 + a2_6*x2*x6 + a2_7*x2*x7 + a2_8*x2*x8 + a2_9*x2*x9 + a2_10*x2*x10 + a2_11*x2*x11 + a2_12*x2*x12;
        dydt(3) = r3*x3 + a3_3*x3^2 + a3_1*x1*x3 + a3_2*x2*x3 + a3_4*x3*x4 + a3_5*x3*x5 + a3_6*x3*x6 + a3_7*x3*x7 + a3_8*x3*x8 + a3_9*x3*x9 + a3_10*x3*x10 + a3_11*x3*x11 + a3_12*x3*x12;
        dydt(4) = r4*x4 + a4_4*x4^2 + a4_1*x1*x4 + a4_2*x2*x4 + a4_3*x3*x4 + a4_5*x4*x5 + a4_6*x4*x6 + a4_7*x4*x7 + a4_8*x4*x8 + a4_9*x4*x9 + a4_10*x4*x10 + a4_11*x4*x11 + a4_12*x4*x12;
        dydt(5) = r5*x5 + a5_5*x5^2 + a5_1*x1*x5 + a5_2*x2*x5 + a5_3*x3*x5 + a5_4*x4*x5 + a5_6*x5*x6 + a5_7*x5*x7 + a5_8*x5*x8 + a5_9*x5*x9 + a5_10*x5*x10 + a5_11*x5*x11 + a5_12*x5*x12;
        dydt(6) = r6*x6 + a6_6*x6^2 + a6_1*x1*x6 + a6_2*x2*x6 + a6_3*x3*x6 + a6_4*x4*x6 + a6_5*x5*x6 + a6_7*x6*x7 + a6_8*x6*x8 + a6_9*x6*x9 + a6_10*x6*x10 + a6_11*x6*x11 + a6_12*x6*x12;
        dydt(7) = r7*x7 + a7_7*x7^2 + a7_1*x1*x7 + a7_2*x2*x7 + a7_3*x3*x7 + a7_4*x4*x7 + a7_5*x5*x7 + a7_6*x6*x7 + a7_8*x7*x8 + a7_9*x7*x9 + a7_10*x7*x10 + a7_11*x7*x11 + a7_12*x7*x12;
        dydt(8) = r8*x8 + a8_8*x8^2 + a8_1*x1*x8 + a8_2*x2*x8 + a8_3*x3*x8 + a8_4*x4*x8 + a8_5*x5*x8 + a8_6*x6*x8 + a8_7*x7*x8 + a8_9*x8*x9 + a8_10*x8*x10 + a8_11*x8*x11 + a8_12*x8*x12;
        dydt(9) = r9*x9 + a9_9*x9^2 + a9_1*x1*x9 + a9_2*x2*x9 + a9_3*x3*x9 + a9_4*x4*x9 + a9_5*x5*x9 + a9_6*x6*x9 + a9_7*x7*x9 + a9_8*x8*x9 + a9_10*x9*x10 + a9_11*x9*x11 + a9_12*x9*x12;
        dydt(10) = r10*x10 + a10_10*x10^2 + a10_1*x1*x10 + a10_2*x2*x10 + a10_3*x3*x10 + a10_4*x4*x10 + a10_5*x5*x10 + a10_6*x6*x10 + a10_7*x7*x10 + a10_8*x8*x10 + a10_9*x9*x10 + a10_11*x10*x11 + a10_12*x10*x12;
        dydt(11) = r11*x11 + a11_11*x11^2 + a11_1*x1*x11 + a11_2*x2*x11 + a11_3*x3*x11 + a11_4*x4*x11 + a11_5*x5*x11 + a11_6*x6*x11 + a11_7*x7*x11 + a11_8*x8*x11 + a11_9*x9*x11 + a11_10*x10*x11 + a11_12*x11*x12;
        dydt(12) = r12*x12 + a12_12*x12^2 + a12_1*x1*x12 + a12_2*x2*x12 + a12_3*x3*x12 + a12_4*x4*x12 + a12_5*x5*x12 + a12_6*x6*x12 + a12_7*x7*x12 + a12_8*x8*x12 + a12_9*x9*x12 + a12_10*x10*x12 + a12_11*x11*x12;

        output = dydt';
    end
end