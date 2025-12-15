function [baseMVA,bus,gen,branch,area,gencost,nomnod,linee,nomare,lim_area,cpiu,cmeno,capa] = Grid_WSCC_Sauer
%CASE14    Power flow data for WSCC Grid, pag. 171 Sauer
%%-----  Power Flow Data  -----%%
%% system MVA base
baseMVA = 100;
%% bus data
% bus_i	type Pd     Qd	    Gs	Bs   area    Vm	    Va    baseKV zone Vmax	Vmin
bus = [
	1	3	0        0       0	0       1	1.04   	0	    110	 1	1.06	0.94;
	2	2	0        0       0	0       1	1.025	0   	110	 1	1.06	0.94;
	3	2	0        0       0	0    	1	1.025	0   	110	 1	1.06	0.94;
    4	1	0        0       0	0    	1	1		0   	110	 1	1.06	0.94;
    5	1	125      50      0	0    	1	1		0   	110	 1	1.06	0.94;
    6	1	90       30      0	0    	1	1		0   	110	 1	1.06	0.94;
    7	1	0        0       0	0    	1	1		0   	110	 1	1.06	0.94;
    8	1	100      35      0	0    	1	1		0   	110	 1	1.06	0.94;
    9	1	0        0       0	0    	1	1		0   	110	 1	1.06	0.94;
];
% % PQ		= 1;
% % PV		= 2;
% % REF		= 3;
% % NONE	= 4;
%% generator data
%	bus	Pg	     Qg	Qmax     Qmin    Vg   mBase	   status	Pmax	Pmin
gen = [
	1	71.6   	 0 	1000	-1000	 1.04	100	       1	90      0;
	2	163	     0 	1000	-1000	 1.025	100	       1	90	    0;
    3	85	     0 	1000	-1000	 1.025	100	       1	90	    0;
];

%% branch data
%	  fbus tbus     r       x           b      rateA	rateB rateC	ratio	angle	status
branch = [
        1   4		0       0.0576    0.0       9900	0       0     0       0     1;
        2   7		0       0.0625	  0.0       9900	0       0	  0       0     1;
        3   9		0       0.0586	  0.0       9900	0       0	  0       0     1;
        4   5		0.010	0.085	  0.088*2	9900	0       0	  0       0     1;
        4   6		0.017   0.092	  0.079*2   9900	0       0	  0       0     1;
        5   7		0.032   0.161	  0.153*2   9900	0       0	  0       0     1;
        6   9		0.039   0.17	  0.179*2   9900	0       0	  0       0     1;
        7   8		0.0085  0.072	  0.0745*2  9900	0       0	  0       0     1;
        8   9		0.0119  0.1008	  0.1045*2  9900	0       0	  0       0     1;
];

%% Others

% Bus names:
nomnod = [
'bus      1'   
'bus      2'   
'bus      3'   
'bus      4'   
'bus      5'   
'bus      6'   
'bus      7'   
'bus      8'   
'bus      9'
];   

% Lines name
linee = ['from   bus      1     to   bus      4',... 

'from   bus      2     to   bus      7',... 

'from   bus      3     to   bus      9',... 

'from   bus      4     to   bus      5',... 

'from   bus      4     to   bus      6',... 

'from   bus      5     to   bus      7',... 

'from   bus      6     to   bus      9',... 

'from   bus      8     to   bus      9',... 
];

%% -----  OPF Data  -----%% - Not used
%% area data
area = [
	1	1;
];

nomare = str2mat('Area1');


%% generator cost data
%	1	startup	shutdown	n	x0	y0	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
gencost = [
	2	0	0	3	0.0430293	20	0;
	2	0	0	3	0.25	20	0;
	2	0	0	3	0.01	40	0;
];

%% active power limit per area

lim_area =[   
   1 220 100
];


cpiu = [];
cmeno = [];
capa =  [];
return;
