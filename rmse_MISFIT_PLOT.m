%% close all mamrousi model
x=1:1:40;
load  SGRDL_M.mat
y1 = out1.R;
y2 = out2.R;
y3 =out3.R;
y4 = out4.R;
y_SGRDL=[y1' y2' y3' y4'];
clear y1 y2 y3 y4
load  ASRI_M
y11 = out1.R;
y12 = out2.R;
y13 = out3.R;
y14 = out4.R;
y_ASRI=[y11' y12' y13' y14'];
load 'TV_M' 
load 'WR_M'
y_TV = outTV.R;
y_WR=outWR.R;
plot(x,y_WR,'k',x,y_TV,'b',x,y_ASRI,'g',x,y_SGRDL,'r')
%% close all
x=1:1:40;
load  SGRDL_M.mat
y1 = out1.JJ;
y2 = out2.JJ;
y3 =out3.JJ;
y4 = out4.JJ;
y_SGRDL=[y1' y2' y3' y4'];
clear y1 y2 y3 y4
load  ASRI_M
y11 = out1.J;
y12 = out2.J;
y13 = out3.J;
y14 = out4.J;
y_ASRI=[y11' y12' y13' y14'];
load 'TV_M' 
load 'WR_M'
y_TV = outTV.JJ;
y_WR=outWR.J;
plot(x,y_WR,'k',x,y_TV,'b',x,y_ASRI,'g',x,y_SGRDL,'r')
y=[y_WR;y_TV;y_ASRI';y_SGRDL'];
[y,ps] = mapminmax(y',0,1);
y_WR=y(1,1:40);
y_TV=y(1,41:80);
y_ASRI=y(1,81:120);
y_SGRDL=y(1,121:160);
plot(x,y_WR,'k',x,y_TV,'b',x,y_ASRI,'g',x,y_SGRDL,'r')

%% close all bg model
x=1:1:40;
load  SGRDL_B.mat
y1 = out1.JJ;
y2 = out2.JJ;
y3 =out3.JJ;
y4 = out4.JJ;
y_SGRDL=[y1' y2' y3' y4'];
clear y1 y2 y3 y4
load  ASRI_B
y11 = out1.J;
y12 = out2.J;
y13 = out3.J;
y14 = out4.J;
y_ASRI=[y11' y12' y13' y14'];
load 'TV_B' 
load 'WR_B'
y_TV = outTV.JJ;
y_WR=outWR.J;
plot(x,y_WR,'k',x,y_TV,'b',x,y_ASRI,'g',x,y_SGRDL,'r')
y=[y_WR;y_TV;y_ASRI';y_SGRDL'];
[y,ps] = mapminmax(y',0,1);
y_WR=y(1,1:40);
y_TV=y(1,41:80);
y_ASRI=y(1,81:120);
y_SGRDL=y(1,121:160);
plot(x,y_WR,'k',x,y_TV,'b',x,y_ASRI,'g',x,y_SGRDL,'r')

%%
x=1:1:40;
load  SGRDL_B.mat
y1 = out1.R;
y2 = out2.R;
y3 =out3.R;
y4 = out4.R;
y_SGRDL=[y1' y2' y3' y4'];
clear y1 y2 y3 y4
load  ASRI_B
y11 = out1.R;
y12 = out2.R;
y13 = out3.R;
y14 = out4.R;
y_ASRI=[y11' y12' y13' y14'];
load 'TV_B' 
load 'WR_B'
y_TV = outTV.R;
y_WR=outWR.R;
plot(x,y_WR,'k',x,y_TV,'b',x,y_ASRI,'g',x,y_SGRDL,'r')
