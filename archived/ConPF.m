%% Continuous Power Flow
function Results=ConPF
clear all; clc;
global bus_data bustype nbus save Result_delta Result_V Result_Sigma Result_lambda Result_iter 

Y = tap_ybus;
Sbase = 100;

Pld = bus_data(:, 7);
Qld = bus_data(:, 8); 

Pgen = bus_data(:, 9);
Qgen = bus_data(:, 10); 

Pinj = (Pgen - Pld)/Sbase;
Qinj = (Qgen - Qld)/Sbase;

%% Initial values for flat start
delta = zeros(nbus,1);
V = ones(nbus,1);

for i=1:nbus
    if bustype(i) == 2 || bustype(i) == 3 
        V(i,1) = bus_data(i, 5);
    end
end

P = zeros(nbus,1);
Q = zeros(nbus,1);
for i=1:nbus
    for j=1:nbus
        P(i) = P(i) + V(i) * V(j) * abs(Y(i,j)) * cos(angle(Y(i,j)) + delta(j) - delta(i));
        Q(i) = Q(i) - V(i) * V(j) * abs(Y(i,j)) * sin(angle(Y(i,j)) + delta(j) - delta(i));
    end
end
% ------------------------------------------------------
%% CPF parameters
save = 0; % This is introduced for saving all results
lambda = 0; 
Section = 1; % Section indicates that we are in increasing p, decreasing V or decreasing P stage of CPF
CPF_bus = 10; % This indicates that on which bus of the power system we are studying the CPF

Pinj_initial = Pinj;
Qinj_initial = Qinj;

q = 0;
for i=2:nbus
    if bustype(i) == 0
       q = q + 1;
       Qinj_initial_(q) = Qinj_initial(i);
    if i == CPF_bus
       CPF_position = q;
    end
    end
end


K = [Pinj_initial(2:end, :); Qinj_initial_']; 

%% CPF code
while lambda >= 0
    J_LF = tap_Jacob(V,delta); 
    CPT_index = nbus - 1 + CPF_position;
    lambda0 = lambda;
    delta0 = delta(2:end);

    q = 0;
    for i=2:nbus
        if bustype(i) == 0
            q = q + 1;
            V0(q) = V(i);
       end
    end

    %% Predictor stage
    J_CPF = [-J_LF, K];
    if Section == 1 
        sigma = 0.1;
        ek = [zeros(1,length(K)),1];
    elseif Section == 2 
        sigma = 0.025;
        ek = [zeros(1,length(K)+1)];
        ek(CPT_index) = -1;
    elseif Section == 3 
        sigma = 0.1;
        ek = [zeros(1,length(K)),-1];
    end
    J_CPF = [J_CPF;ek];
    Predicted_Vector = [delta0; V0'; lambda0] + sigma * inv(J_CPF) * [zeros(length(K),1);1];
    lambda = Predicted_Vector(end);
    V_CPF_bus_Initial = Predicted_Vector(CPT_index);

    %% Corrector Section
    eps = 0.01;
    err = 10;
    iter = 0;
    V_Initial = V;
    delta_Initial = delta;
    
    while err > eps
        P = zeros(nbus,1);
        Q = zeros(nbus,1);
        for i=1:nbus
            for j=1:nbus
                P(i) = P(i) + V(i) * V(j) * abs(Y(i,j)) * cos(angle(Y(i,j)) + delta(j) - delta(i));
                Q(i) = Q(i) - V(i) * V(j) * abs(Y(i,j)) * sin(angle(Y(i,j)) + delta(j) - delta(i));
            end
        end
        J_LF = tap_Jacob(V,delta);
        Pinj = lambda.*Pinj_initial;
        Qinj = lambda.*Qinj_initial;
        deltaP = Pinj(2:nbus) - P(2:nbus);

        k = 0;
        for i=2:nbus
            if bustype(i) == 0
                k = k + 1;
                deltaQ(k) = Qinj(i) - Q(i);
            end
        end
        if Section == 1 || Section == 3
            deltaP_Q = cat(1, deltaP, deltaQ');
            MisMatch = SolveWithLU(J_LF,deltaP_Q);
        elseif Section == 2
            del_V_new = V_CPF_bus_Initial - V(CPF_bus);
            deltaP_Q_V = [deltaP; deltaQ'; del_V_new];
            J_LF_New = [-J_LF,K];
            J_LF_New = [J_LF_New; ek];
            MisMatch = -SolveWithLU(J_LF_New,deltaP_Q_V);
            lambda = lambda + MisMatch(end);
        end
        del_delta = MisMatch(1:nbus-1);
        del_V = MisMatch(nbus:end);
        delta(2:end) = delta(2:end) + del_delta;

        k = 0;
        for i=2:nbus
            if bustype(i) == 0
                k = k + 1;
                V(i) = V(i) + V(i)*del_V(k);
            end
        end

        %% Calculating new set of P and Q with updated V and deltata
        P = zeros(nbus,1);
        Q = zeros(nbus,1);
        for i=1:nbus
            for j=1:nbus
                P(i) = P(i) + V(i) * V(j) * abs(Y(i,j)) * cos(angle(Y(i,j)) + delta(j) - delta(i));
                Q(i) = Q(i) - V(i) * V(j) * abs(Y(i,j)) * sin(angle(Y(i,j)) + delta(j) - delta(i));
            end
        end
        deltaP = Pinj(2:nbus) - P(2:nbus);

        k=0;
        for i=2:nbus
            if bustype(i) == 0
                k = k + 1;
                deltaQ(k) = Qinj(i)-Q(i);
            end
        end

        err = max(abs(cat(1, deltaP, deltaQ')));
        if iter > 6
            V = V_Initial;
            delta = delta_Initial;
            lambda = lambda0;
            if Section == 1
                Section = 2;
                break
            elseif Section == 2
                Section =3;
                break
            end
        end
        iter = iter + 1;
    end
    save = save + 1;
    Result_iter(save) = iter;
    Result_lambda(save) = lambda;
    Result_Sigma(save) = sigma;
    Result_V(:,save) = V;
    Result_delta(:,save) = delta;

end
clc
CPF_Plot;
disp('CPF increments:');
disp('------------------');
disp('  Number | lambda | Bus.10 Vol | Sigma');
Results=[(1:save-1)',Result_lambda(1:save-1)',Result_V(9,(1:save-1))',Result_Sigma(1:save-1)'];
end
function CPF_Plot
global nbus  Result_V  Result_lambda bustype
for i=2:nbus
        if bustype(i) == 0
            figure(i)
            plot(Result_lambda,Result_V(i,:),'*');
            xlim([0 4.5]);
            hold on;
            xlabel('lambda (Change in power)');
            ylabel(['Bus.',num2str(i)]);
            title(['CPF for bus',num2str(i)]);
        end
   
end

end
