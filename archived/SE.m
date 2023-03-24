function V_delta = SE

clc; clear all;

global bus_data nbus line_data V_meas delta_meas V_var PQ_var Hx Pinj_meas ...
       Qinj_meas Pij_meas Pji_meas Qij_meas Qji_meas;

Y = tap_ybus;
Sbase = 100;

Pld = bus_data(:, 7);
Qld = bus_data(:, 8); 

Pgen = bus_data(:, 9);
Qgen = bus_data(:, 10); 

Pinj = (Pgen - Pld)/Sbase;
Qinj = (Qgen - Qld)/Sbase;

% Initial values for flat start

delta = bus_data(:,6);  % These are in degrees
V(:, 1) = bus_data(:, 5);

fb = line_data(:,1);     % From bus number
tb = line_data(:,2);     % To bus number

[d1,~] = size(line_data);                  
nbranch = d1;            % Number of branches
b = line_data(:,9);      % Ground Admittance, B/2

%% initiating the measurements
V_meas = V;
delta_meas = delta;
Pinj_meas = Pinj;
Qinj_meas = Qinj;

for i = 1:nbranch
    Iij_meas(i, 1) = (V_meas(fb(i)) - V_meas(tb(i))) * Y(fb(i), tb(i));

    Pij_meas(i,1) = -real(Y(fb(i), tb(i))) * (V_meas(fb(i))^2) + ...
        abs(V_meas(fb(i)) * V_meas(tb(i)) * Y(fb(i), tb(i))) * cos(angle(Y(fb(i), tb(i))) + ...
        delta_meas(tb(i)) - delta_meas(fb(i)));

    Qij_meas(i,1) = -(-imag(Y(fb(i), tb(i))) + abs(b(tb(i)))/2) * abs(V_meas(fb(i))^2) ...
        - abs(V_meas(fb(i)) * V_meas(tb(i)) * Y(fb(i),tb(i))) * sin(angle(Y(fb(i), tb(i))) ... 
        + delta_meas(tb(i)) - delta_meas(fb(i)));

    Pji_meas(i,1) = -real(Y(tb(i),fb(i))) * (V_meas(tb(i))^2) + ...
        abs(V_meas(tb(i)) * V_meas(fb(i)) * Y(tb(i), fb(i))) * cos(angle(Y(tb(i), fb(i))) ...
        + delta_meas(fb(i)) - delta_meas(tb(i)));

    Qji_meas(i,1) = -(-imag(Y(tb(i), fb(i)))+abs(b(fb(i)))/2) * abs(V_meas(tb(i))^2) - ...
        abs(V_meas(tb(i)) * V_meas(fb(i)) * Y(tb(i), fb(i))) * sin(angle(Y(tb(i), fb(i))) ...
        + delta_meas(fb(i)) - delta_meas(tb(i)));
end

% adding noise to the measurements
V_var = 0.02;
PQ_var = 0.01;

V_meas = V_meas + V_meas * randn * V_var;
Pinj_meas = Pinj_meas + Pinj_meas * randn * PQ_var;
Qinj_meas = Qinj_meas + Qinj_meas * randn * PQ_var;
Pij_meas = Pij_meas + Pij_meas * randn * PQ_var;
Pji_meas = Pji_meas + Pji_meas * randn * PQ_var;
Qij_meas = Qij_meas + Qij_meas * randn * PQ_var;
Qji_meas = Qji_meas + Qji_meas * randn * PQ_var;

%% initiating the state estimation

eps = 0.001;
err = 1;
iter = 0;

% introducing some bad data 
V_meas(5) = 4;
Pij_meas(5) = 20;
Qij_meas(9) = 20;



    while err >= eps
        % Pbus and Qbus
        P = zeros(nbus, 1);
        Q = zeros(nbus, 1);
        for i=1:nbus
            for j=1:nbus
                P(i) = P(i) + V(i) * V(j) * abs(Y(i,j)) * cos(angle(Y(i,j)) + delta(j) - delta(i));
                Q(i) = Q(i) - V(i) * V(j) * abs(Y(i,j)) * sin(angle(Y(i,j)) + delta(j) - delta(i));
            end
        end

        % Pline and Qline
        for i = 1:nbranch
                Pij(i,1) = -real(Y(fb(i), tb(i))) * abs(V(fb(i))^2) ...
                    + abs(V(fb(i)) * V(tb(i)) * Y(fb(i), tb(i))) * cos(angle(Y(fb(i), tb(i))) ...
                    + delta(tb(i)) - delta(fb(i)));

                Qij(i,1) = -(-imag(Y(fb(i), tb(i))) + abs(b(tb(i)))/2) * abs(V(fb(i))^2) ...
                    - abs(V(fb(i)) * V(tb(i)) * Y(fb(i),tb(i))) * sin(angle(Y(fb(i), tb(i))) ...
                    + delta(tb(i)) - delta(fb(i)));

                Pji(i,1) = -real(Y(tb(i), fb(i))) * abs(V(tb(i))^2) ...
                    + abs(V(tb(i)) * V(fb(i)) * Y(tb(i), fb(i))) * cos(angle(Y(tb(i), fb(i))) ...
                    + delta(fb(i)) - delta(tb(i)));

                Qji(i,1) = -(-imag(Y(tb(i), fb(i))) + abs(b(fb(i)))/2) * abs(V(tb(i))^2) ...
                    - abs(V(tb(i)) * V(fb(i)) * Y(tb(i), fb(i))) * sin(angle(Y(tb(i), fb(i))) ...
                    + delta(fb(i)) - delta(tb(i)));
        end

        % forming the error vector
        Err_vec = [V_meas - V ; Pinj_meas - Pinj ; Qinj_meas - Qinj ; Pij_meas - Pij ;...
                   Pji_meas - Pji ; Qij_meas - Qij ; Qji_meas - Qji];
        % Forming the matrix of partial derivatives
        Hx = H_x(V, delta, P, Q);

        % Weight matrix
        R = (diag([ones(1, length(V)) * V_var, ones(1,length(Hx) - length(V)) * PQ_var]))^2;

        % minimizing the error
        Xhat = [delta(2:end); V]; 
        Xhat_Previous = Xhat;
        Xhat = Xhat + inv(Hx' * inv(R) * Hx) * Hx' * inv(R) * Err_vec;
        delta(2:end) = Xhat(1:nbus-1);
        V = Xhat(nbus:end);
        iter = iter + 1;
        err = max(abs(Xhat - Xhat_Previous));
    end

disp('-----------------------');
disp('    Voltages & angles:');
V_delta = [V, delta];

P = zeros(nbus,1);
Q = zeros(nbus,1);
for i=1:nbus
    for j=1:nbus
        P(i) = P(i) + V(i) * V(j) * abs(Y(i,j)) * cos(angle(Y(i,j)) + delta(j) - delta(i));
        Q(i) = Q(i) - V(i) * V(j) * abs(Y(i,j)) * sin(angle(Y(i,j)) + delta(j) - delta(i));
    end
end

for i = 1:nbranch
    Pij(i,1) = -real(Y(fb(i),tb(i)))*(V_meas(fb(i))^2) + abs(V_meas(fb(i))*V_meas(tb(i))*Y(fb(i),tb(i)))*cos(angle(Y(fb(i),tb(i)))+delta_meas(tb(i))-delta_meas(fb(i)));
    Qij(i,1) = -(-imag(Y(fb(i),tb(i)))+abs(b(tb(i)))/2)*abs(V_meas(fb(i))^2) - abs(V_meas(fb(i))*V_meas(tb(i))*Y(fb(i),tb(i)))*sin(angle(Y(fb(i),tb(i)))+delta_meas(tb(i))-delta_meas(fb(i)));
    Pji(i,1) = -real(Y(tb(i),fb(i)))*(V_meas(tb(i))^2) + abs(V_meas(tb(i))*V_meas(fb(i))*Y(tb(i),fb(i)))*cos(angle(Y(tb(i),fb(i)))+delta_meas(fb(i))-delta_meas(tb(i)));
    Qji(i,1) = -(-imag(Y(tb(i),fb(i)))+abs(b(fb(i)))/2)*abs(V_meas(tb(i))^2) - abs(V_meas(tb(i))*V_meas(fb(i))*Y(tb(i),fb(i)))*sin(angle(Y(tb(i),fb(i)))+delta_meas(fb(i))-delta_meas(tb(i)));
end

Err_vec = [V_meas - V ; Pinj_meas - Pinj ; Qinj_meas - Qinj ; Pij_meas - Pij ;...
           Pji_meas - Pji ; Qij_meas - Qij ; Qji_meas - Qji];

f = 0;
for i = 1:length(Err_vec)
    f = f + R(i, i) * Err_vec(i).^2;
end

disp('-----------------------');
disp('Error is:');
disp(f)

k = length(Err_vec) - length(Xhat);
Err_thresh = chi2inv(0.99, k)

%if f >= Err_thresh
%    disp('Bad data detected')
%else
%    disp('Measurements are good')

errors = (V - bus_data(:, 5))./bus_data(:, 5) * 100;
%plot(errors)

end