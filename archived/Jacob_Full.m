%% Calculation Jacobian Matrix
function J = Jacob_Full(V, delta, P, Q)
global nbus 

Y = tap_ybus;

for m = 1:nbus 
    for n = 1:nbus
        if m ~= n
            J11(m, n) = -abs(V(m) * V(n) * Y(m, n)) * sin(angle(Y(m, n)) + delta(n) - delta(m));
            J21(m, n) = -abs(V(m) * V(n) * Y(m, n)) * cos(angle(Y(m, n)) + delta(n) - delta(m));
            J12(m, n) = abs(V(m) * V(n) * Y(m, n)) * cos(angle(Y(m, n)) + delta(n) - delta(m));
            J22(m, n) = -abs(V(m) * V(n) * Y(m, n)) * sin(angle(Y(m, n)) + delta(n) - delta(m));
        else
            J11(m, m) = -Q(m) - (abs(V(m))^2) * imag(Y(m,m));
            J21(m, m) = P(m) - (abs(V(m))^2) * real(Y(m,m));
            J12(m, m) = P(m) + (abs(V(m))^2) * real(Y(m,m));
            J22(m, m) = Q(m) - (abs(V(m))^2) * imag(Y(m,m));
        end
    end
end
J=[J11 J12; J21 J22];

end