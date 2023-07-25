function [J, JTable] = constructJacobian(P, Q, ...
    V, delta, N, ybus, BMatrix, ~, ...
    nPQ, nPV, listOfPQBuses, listOfNonSlackBuses, powerFlowMethod, desiredOutput)
    if strcmp(powerFlowMethod, 'NRPF') || strcmp(powerFlowMethod, 'Decoupled NRPF')
        if strcmp(desiredOutput, 'both')
            J11 = zeros(nPQ+nPV, nPQ+nPV);
            for row = 1:nPQ+nPV
                for col = 1:nPQ+nPV
                    i = listOfNonSlackBuses(row);
                    k = listOfNonSlackBuses(col);
                    if k == i
                        J11(row, col) = -Q(i) - imag(ybus(i,i))*abs(V(i))^2;
                    else
                        J11(row, col) = -abs(ybus(i,k))*abs(V(k))*abs(V(i))*sin( angle(ybus(i,k)) + delta(k) - delta(i) );
                    end
                end
            end
        
            J12 = zeros(nPQ+nPV, nPQ);
            for row = 1:nPQ+nPV
                for col = 1:nPQ
                    i = listOfNonSlackBuses(row);
                    k = listOfPQBuses(col);
                    if k == i
                        J12(row, col) = P(i) + real(ybus(i, i))*abs(V(i))^2;
                    else
                        J12(row, col) = abs(ybus(i, k))*abs(V(k))*abs(V(i))*cos( angle(ybus(i,k)) + delta(k) - delta(i) );
                    end
                end
            end
            
            J21 = zeros(nPQ, nPQ+nPV);
            for row = 1:nPQ
                for col = 1:nPQ+nPV
                    i = listOfPQBuses(row);
                    k = listOfNonSlackBuses(col);
                    if k == i
                        J21(row, col) = P(i) - real(ybus(i, i))*abs(V(i))^2;
                    else
                        J21(row, col) = -abs(ybus(i, k))*abs(V(k))*abs(V(i))*cos( angle( ybus(i, k) ) + delta(k) - delta(i) );
                    end
                end
            end
        
            J22 = zeros(nPQ, nPQ);
            for row = 1:nPQ
                for col = 1:nPQ
                    i = listOfPQBuses(row);
                    k = listOfPQBuses(col);
                    if k == i
                        J22(row, col) = Q(i) - imag(ybus(i, i))*abs(V(i))^2;
                    else
                        J22(row, col) = -abs(ybus(i, k))*abs(V(k))*abs(V(i))*sin( angle( ybus(i, k)) + delta(k) - delta(i));
                    end
                end
            end
        
            J = [J11 J12; J21 J22];
        
            jacobianNamesP = cell(N, 1);
            jacobianNamesQ = cell(N, 1);
            jacobianNamesDelta = cell(N, 1);
            jacobianNamesDeltaVByV = cell(N, 1);
            for row = 1:N
                jacobianNamesP(row) = {sprintf('$P_%i$', row)};
                jacobianNamesDelta(row) = {sprintf('$delta_%i$', row)};
                jacobianNamesQ(row) = {sprintf('$Q_%i$', row)};
                jacobianNamesDeltaVByV(row) = {sprintf('$DeltaVByV_%i$', row)};
            end
            jacobianNamesP = jacobianNamesP(listOfNonSlackBuses);
            jacobianNamesDelta = jacobianNamesDelta(listOfNonSlackBuses);
            jacobianNamesQ = jacobianNamesQ(listOfPQBuses);
            jacobianNamesDeltaVByV = jacobianNamesDeltaVByV(listOfPQBuses);
        
            JTable = array2table(J, ...
                RowNames=[jacobianNamesP; jacobianNamesQ],...
                VariableNames = [jacobianNamesDelta; jacobianNamesDeltaVByV]);
        elseif strcmp(desiredOutput, 'P')
            J11 = zeros(nPQ+nPV, nPQ+nPV);
            for row = 1:nPQ+nPV
                for col = 1:nPQ+nPV
                    i = listOfNonSlackBuses(row);
                    k = listOfNonSlackBuses(col);
                    if k == i
                        J11(row, col) = -Q(i) - imag(ybus(i,i))*abs(V(i))^2;
                    else
                        J11(row, col) = -abs(ybus(i,k))*abs(V(k))*abs(V(i))*sin( angle(ybus(i,k)) + delta(k) - delta(i) );
                    end
                end
            end
            J = J11;
            jacobianNamesP = cell(N, 1);
            jacobianNamesDelta = cell(N, 1);
            for row = 1:N
                jacobianNamesP(row) = {sprintf('$P_%i$', row)};
                jacobianNamesDelta(row) = {sprintf('$delta_%i$', row)};
            end
            jacobianNamesP = jacobianNamesP(listOfNonSlackBuses);
            jacobianNamesDelta = jacobianNamesDelta(listOfNonSlackBuses);
        
            JTable = array2table(J, ...
                RowNames=jacobianNamesP,...
                VariableNames = jacobianNamesDelta);
        
        elseif strcmp(desiredOutput, 'Q')
            J22 = zeros(nPQ, nPQ);
            for row = 1:nPQ
                for col = 1:nPQ
                    i = listOfPQBuses(row);
                    k = listOfPQBuses(col);
                    if k == i
                        J22(row, col) = Q(i) - imag(ybus(i, i))*abs(V(i))^2;
                    else
                        J22(row, col) = -abs(ybus(i, k))*abs(V(k))*abs(V(i))*sin( angle( ybus(i, k)) + delta(k) - delta(i));
                    end
                end
            end
        
            J = J22;
        
            jacobianNamesQ = cell(N, 1);
            jacobianNamesDeltaVByV = cell(N, 1);
            for row = 1:N
                jacobianNamesQ(row) = {sprintf('$Q_%i$', row)};
                jacobianNamesDeltaVByV(row) = {sprintf('$DeltaVByV_%i$', row)};
            end
            jacobianNamesQ = jacobianNamesQ(listOfPQBuses);
            jacobianNamesDeltaVByV = jacobianNamesDeltaVByV(listOfPQBuses);
        
            JTable = array2table(J, ...
                RowNames=jacobianNamesQ,...
                VariableNames = jacobianNamesDeltaVByV);
        end
    elseif strcmp(powerFlowMethod, 'Fast Decoupled NRPF')
        ybus = BMatrix;
        if strcmp(desiredOutput, 'P')
            J11 = zeros(nPQ+nPV, nPQ+nPV);
            for row = 1:nPQ+nPV
                for col = 1:nPQ+nPV
                    i = listOfNonSlackBuses(row);
                    k = listOfNonSlackBuses(col);
                    if k == i
                        J11(row, col) = -Q(i) - ybus(i,i)*abs(V(i))^2;
                    else
                        J11(row, col) = -ybus(i,k)*abs(V(k))*abs(V(i));
                    end
                end
            end
            J = J11;
            jacobianNamesP = cell(N, 1);
            jacobianNamesDelta = cell(N, 1);
            for row = 1:N
                jacobianNamesP(row) = {sprintf('$P_%i$', row)};
                jacobianNamesDelta(row) = {sprintf('$delta_%i$', row)};
            end
            jacobianNamesP = jacobianNamesP(listOfNonSlackBuses);
            jacobianNamesDelta = jacobianNamesDelta(listOfNonSlackBuses);
        
            JTable = array2table(J, ...
                RowNames=jacobianNamesP,...
                VariableNames = jacobianNamesDelta);
        
        elseif strcmp(desiredOutput, 'Q')
            J22 = zeros(nPQ, nPQ);
            for row = 1:nPQ
                for col = 1:nPQ
                    i = listOfPQBuses(row);
                    k = listOfPQBuses(col);
                    if k == i
                        J22(row, col) = Q(i) - ybus(i, i)*abs(V(i))^2;
                    else
                        J22(row, col) = -ybus(i, k)*abs(V(k))*abs(V(i))*( delta(k) - delta(i));
                    end
                end
            end
        
            J = J22;
        
            jacobianNamesQ = cell(N, 1);
            jacobianNamesDeltaVByV = cell(N, 1);
            for row = 1:N
                jacobianNamesQ(row) = {sprintf('$Q_%i$', row)};
                jacobianNamesDeltaVByV(row) = {sprintf('$DeltaVByV_%i$', row)};
            end
            jacobianNamesQ = jacobianNamesQ(listOfPQBuses);
            jacobianNamesDeltaVByV = jacobianNamesDeltaVByV(listOfPQBuses);
        
            JTable = array2table(J, ...
                RowNames=jacobianNamesQ,...
                VariableNames = jacobianNamesDeltaVByV);
        end
    end
end