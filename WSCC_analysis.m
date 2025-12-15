clear
close all
clc
warning off
path(path,'Load_flow');
%                                                                         AVR   PSS  
[A, B, pos_var, ~] = small_sgn_linp('Grid_WSCC_Sauer','gendat_WSCC_9_bus',  1,  0);
