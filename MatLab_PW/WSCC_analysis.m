clear
close all
clc
warning off
path(path,'Load_flow');
%                                                                         AVR   PSS  
[A, B, pos_var, ~] = small_sgn_linp('Grid_WSCC_Sauer','gendat_WSCC_9_bus',  0,  0);


% Perform eigenvalue analysis on the system matrix A
[V, eigenValues] = eig(A);
W = inv(V).';
eigenValues = diag(eigenValues);
realIndex = abs(imag(eigenValues)) < 5e-3 ; %vettori logico che indica gli auitovalori reali

figure
plot(real(eigenValues(realIndex)), imag(eigenValues(realIndex)), 'o', 'MarkerSize',8, 'LineWidth',1.5);
hold on;
plot(real(eigenValues(~realIndex)), imag(eigenValues(~realIndex)), 'x', 'MarkerSize',8, 'LineWidth',1.5)
axis equal
grid on
xlabel('Re(z)')
ylabel('Im(z)')
title('Autovalori nel piano complesso')

V_useful = V([4 5 6], ~realIndex); %prendo solo quelli corrisppondenti a autovalori complessi
V_useful = V_useful(:, 1:2:end); %prendo solo uno dei due complessi coniugati
[maxVal, idx] = max(abs(V_useful));


figure
C1 = compassplot(((V_useful(:, 1))./maxVal(1)*exp(-1j*angle(V_useful(idx(1),1)))));

figure
C3 = compassplot(((V_useful(:, 2))./maxVal(2)*exp(-1j*angle(V_useful(idx(2), 2)))));