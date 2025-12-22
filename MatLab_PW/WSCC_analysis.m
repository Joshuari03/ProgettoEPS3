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
realIndex = abs(imag(eigenValues)) < 5e-3 ; %vettore logico che indica gli autovalori reali

figure
plot(real(eigenValues(realIndex)), imag(eigenValues(realIndex)), 'o', 'MarkerSize',8, 'LineWidth',1.5);
hold on;
plot(real(eigenValues(~realIndex)), imag(eigenValues(~realIndex)), 'x', 'MarkerSize',8, 'LineWidth',1.5)
axis equal
grid on
xlabel('Re(z)')
ylabel('Im(z)')
title('Autovalori nel piano complesso')
hold off

V_useful = V([4 5 6], ~realIndex); %prendo solo gli elementi corrisppondenti a autovalori complessi
                                   %e relativi alle omega dei generatori
V_useful = V_useful(:, 1:2:end); %prendo solo uno dei due complessi coniugati
[maxVal, idx] = max(abs(V_useful));

%% Plot modshapes
C = ((V_useful(:, 1))./maxVal(1)*exp(-1j*angle(V_useful(idx(1),1))));
h = gobjects(numel(C),1);
figure;
ax = polaraxes; hold (ax, 'on');
names = ["G1", "G2", "G3"];
for i = 1:numel(C)
    h(i) = compassplot(C(i));
    h(i).DisplayName = names(i);
end
legend (ax, 'show')
hold (ax, 'off')

C = ((V_useful(:, 2))./maxVal(2)*exp(-1j*angle(V_useful(idx(2), 2))));
figure;
ax2 = polaraxes; hold (ax2, 'on');
h = gobjects(numel(C),1);
for i = 1:numel(C)
    h(i) = compassplot(C(i));
    h(i).DisplayName = names(i);
end
legend (ax2, 'show')
hold (ax2, 'off')
%% damping ratio and oscillation drequency
alpha = real(eigenValues(~realIndex));
alpha = alpha(1:2:end);
omega_d = imag(eigenValues(~realIndex));
omega_d = omega_d(1:2:end);
fd = omega_d/(2*pi);

damping_ratio = cos(pi-angle(eigenValues(~realIndex)));
damping_ratio = damping_ratio(1:2:end);

zeta = -alpha./sqrt((alpha).^2 +omega_d.^2);
%% Part. factors
P = V.*conj(W);