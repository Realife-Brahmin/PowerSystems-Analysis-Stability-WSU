function Q_inj_1_Formula = calculateReactivePowerInjection(x, V, ybus, E)
    Q_inj_1_Formula = -imag(ybus(1, 1))*abs(V(1)^2);
    listOfNeighbours = E{1};
    for k = listOfNeighbours
        if k == 5
            Vk = x(15);
        else
            Vk = V(k);
        end
        Q_inj_1_Formula = Q_inj_1_Formula - ...
            abs( ybus(1, k) * Vk * V(1) )*...
            sin( angle(ybus(1, k)) + x(k-1));
    end
end