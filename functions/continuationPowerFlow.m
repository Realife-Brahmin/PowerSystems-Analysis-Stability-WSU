function [V_CPF, delta_CPF, lambda_CPF, iter_CPF, sigma_CPF] =...
    continuationPowerFlow(busData, ybus, PSpec, QSpec, V, delta, BMatrix, E, ...
    nPQ, nPV, listOfPQBuses, listOfNonSlackBuses, powerFlowMethod, desiredOutput, CPF_Bus, displayCPFResults, plotCPFPlots)
    
    N = size(busData, 1);
    P = zeros(N, 1);
    Q = zeros(N, 1);

    CPFItrMax = 100;
    sectionItrMax = 12;
    verbose = false;
    sigma1 = 0.1;
    sigma2 = 0.005;
    sigma3 = 0.1;

    for i = 1 : N
        for k = 1 : N
            P(i) = P(i) + V(i) * V(k) * abs(ybus(i,k)) * cos( angle( ybus(i,k) ) + delta(k) - delta(i) );
            Q(i) = Q(i) - V(i) * V(k) * abs(ybus(i,k)) * sin( angle( ybus(i,k) ) + delta(k) - delta(i) );
        end
    end

    storeIdx = 0;
    lambda = 0; 
    Section = 1; % Section indicates that we are in increasing p, decreasing V or decreasing P stage of CPF

    P0 = PSpec;
    Q0 = QSpec;

    Q_PQ = QSpec(listOfPQBuses);
    CPF_position = find(listOfPQBuses == CPF_Bus);

    
    K = [P0(listOfNonSlackBuses, :); Q_PQ]; 

    [V_CPF, delta_CPF] = deal(zeros(N, CPFItrMax));
    [iter_CPF, lambda_CPF, sigma_CPF] = deal(zeros(CPFItrMax, 1));

%% Begin Continuation Power Flow Iterations

    while lambda >= 0
% Stage: Initial Power Flow and basic housekeeping

        [J_LF, ~] = constructJacobian(P, Q, ...
            V, delta, N, ybus, BMatrix, E, ...
            nPQ, nPV, listOfPQBuses, listOfNonSlackBuses, powerFlowMethod, desiredOutput);

        CPT_index = N - 1 + CPF_position;
        lambda0 = lambda;
        delta0 = delta(listOfNonSlackBuses);

        V_PQ_0 = V(listOfPQBuses);
% Stage: Predictor 

        J_CPF0 = [-J_LF, K];
        if Section == 1 
            sigma = sigma1;
            ek = [zeros(1,length(K)),1];
        elseif Section == 2 
            sigma = sigma2;
            ek = zeros(1, length(K)+1);
            ek(CPT_index) = -1;
        elseif Section == 3 
            sigma = sigma3;
            ek = [zeros(1, length(K)), -1];
        end

        J_CPF = [J_CPF0; ek];
        b = [zeros(length(K), 1); 1];
        A = J_CPF;
        invAb = A\b; %inv(A)*b
        Predicted_Vector = [delta0; V_PQ_0; lambda0] + sigma * invAb;
        lambda = Predicted_Vector(end);
        V_CPF_bus_Initial = Predicted_Vector(CPT_index);
    
% Stage: Corrector

        eps = 0.01;
        err = 10;
        iter = 0;
        V_Initial = V;
        delta_Initial = delta;
        
        % itr = 1;
        while err > eps

            P = zeros(N, 1);
            Q = zeros(N, 1);

            for i = 1 : N
                for j = 1 : N
                    P(i) = P(i) + V(i) * V(j) * abs(ybus(i,j)) * cos(angle(ybus(i,j)) + delta(j) - delta(i));
                    Q(i) = Q(i) - V(i) * V(j) * abs(ybus(i,j)) * sin(angle(ybus(i,j)) + delta(j) - delta(i));
                end
            end

            [J_LF, ~] = constructJacobian(P, Q, ...
                V, delta, N, ybus, BMatrix, E, ...
                nPQ, nPV, listOfPQBuses, listOfNonSlackBuses, powerFlowMethod, desiredOutput);

            PSpec = lambda .* P0;
            QSpec = lambda .* Q0;
            deltaP = PSpec(listOfNonSlackBuses) - P(listOfNonSlackBuses);
            deltaQ = QSpec(listOfPQBuses) - Q(listOfPQBuses);

            if Section == 1 || Section == 3
                deltaPQ = cat(1, deltaP, deltaQ);
                mismatch = solveUsingLU(J_LF, deltaPQ, iter, verbose);

            elseif Section == 2
                del_V_new = V_CPF_bus_Initial - V(CPF_Bus);
                deltaP_Q_V = [deltaP; deltaQ; del_V_new];
                J_LF_New = [[-J_LF, K]; ek];
                mismatch = -solveUsingLU(J_LF_New, deltaP_Q_V, iter, verbose);
                lambda = lambda + mismatch(end);

            end

            del_delta = mismatch(1 : N-1);
            del_V = mismatch(N : end);
            delta(listOfNonSlackBuses) = delta(listOfNonSlackBuses) + del_delta;
    
            V(listOfPQBuses) = V(listOfPQBuses) .* (ones(nPQ, 1) + del_V(1 : nPQ));

            %% Calculating new set of P and Q with updated V and delta
            P = zeros(N, 1);
            Q = zeros(N, 1);

            for i = 1 : N
                for j = 1 : N
                    P(i) = P(i) + V(i) * V(j) * abs(ybus(i,j)) * cos(angle(ybus(i,j)) + delta(j) - delta(i));
                    Q(i) = Q(i) - V(i) * V(j) * abs(ybus(i,j)) * sin(angle(ybus(i,j)) + delta(j) - delta(i));
                end
            end

            deltaP = PSpec(listOfNonSlackBuses) - P(listOfNonSlackBuses);
            deltaQ = QSpec(listOfPQBuses) - Q(listOfPQBuses);
    
            err = max(abs(cat(1, deltaP, deltaQ)));

            if iter > sectionItrMax
                V = V_Initial;
                delta = delta_Initial;
                lambda = lambda0;
                if Section == 1
                    Section = 2;
                    break
                elseif Section == 2
                    Section = 3;
                    break
                end
           end
            
           iter = iter + 1;
            
        end
        
        storeIdx = storeIdx + 1;
        iter_CPF(storeIdx) = iter;
        lambda_CPF(storeIdx) = lambda;
        sigma_CPF(storeIdx) = sigma;
        V_CPF(:,storeIdx) = V;
        delta_CPF(:, storeIdx) = delta;
    end
    
    iter_CPF = iter_CPF(1:storeIdx);
    lambda_CPF = lambda_CPF(1:storeIdx);
    sigma_CPF = sigma_CPF(1:storeIdx);
    V_CPF = V_CPF(:, 1:storeIdx);
    delta_CPF = delta_CPF(:, 1:storeIdx);

    CPF_Progress = [(1:storeIdx)', lambda_CPF, V_CPF(CPF_Bus, :)', sigma_CPF];
    CPF_Progress_Table = array2table(CPF_Progress, 'VariableNames', {'Itr', 'lambda', strcat('V_', num2str(CPF_Bus)), 'sigma'});

    if plotCPFPlots
        plotCPFResults(V_CPF, lambda_CPF, CPF_Bus, listOfPQBuses);
    end

    if displayCPFResults
        display(CPF_Progress_Table);
    end
end