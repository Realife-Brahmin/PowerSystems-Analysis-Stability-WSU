%% Calculating H_Matrix
function H = H_x(V, delta, P, Q)

global nbus line_data

fb = line_data(:,1);     % From bus number
tb = line_data(:,2);     % To bus number

[d1,~] = size(line_data);                  
nbranch = d1;            % Number of branches
b = line_data(:,9);      % Ground Admittance, B/2

J = Full_Jacob(V, delta, P, Q);
Y = tap_ybus;

%% calculating Hx segments

%% dV/ddelta
H1 = zeros(nbus, nbus-1);

%% dV/dV
H2 = eye(nbus, nbus);

%% dPinj/ddelta
H3 = J(1:nbus, 2:nbus);

%% dPinj/dV
H4 = J(1:nbus, nbus+1:end);

%% dQinj/ddelta
H5 = J(nbus+1:end, 2:nbus);

%% dQinj/dV
H6 = J(nbus+1:end, nbus+1:end);

%% dPij/ddelta
for i = 1:nbranch
    if fb(i)~=1
        H7(i, fb(i)-1) = abs(V(fb(i)) * V(tb(i)) * Y(fb(i), tb(i))) * sin(angle(Y(fb(i), tb(i))) + delta(tb(i)) - delta(fb(i)));
    end
    if   tb(i)~=1
        H7(i, tb(i)-1) = -abs(V(fb(i)) * V(tb(i)) * Y(fb(i),tb(i))) * sin(angle(Y(fb(i), tb(i))) + delta(tb(i)) - delta(fb(i)));
    end
end

%% dPij/dV
for i = 1:nbranch
    H8(i, tb(i)) = abs(V(fb(i)) * Y(fb(i),tb(i))) * cos(angle(Y(fb(i), tb(i))) + delta(tb(i)) - delta(fb(i)));
    H8(i, fb(i)) = -2*V(fb(i)) * real(Y(fb(i), tb(i))) + abs(V(tb(i)) * Y(fb(i), tb(i))) * cos(angle(Y(fb(i), tb(i))) + delta(tb(i)) - delta(fb(i)));
end

%% dPji/ddelta
for i = 1:nbranch
    if tb(i)~=1
        H9(i, tb(i)-1) = abs(V(tb(i)) * V(fb(i)) * Y(tb(i), fb(i))) * sin(angle(Y(tb(i), fb(i))) + delta(fb(i)) - delta(tb(i)));
    end
    if   fb(i)~=1
        H9(i, fb(i)-1) = -abs(V(tb(i)) * V(fb(i)) * Y(tb(i), fb(i))) * sin(angle(Y(tb(i), fb(i))) + delta(fb(i)) - delta(tb(i)));
    end
end

%% dPji/dV
for i = 1:nbranch
    H10(i,fb(i)) = abs(V(tb(i)) * Y(tb(i), fb(i))) * cos(angle(Y(tb(i), fb(i))) + delta(fb(i)) - delta(tb(i)));
    H10(i,tb(i)) = -2*V(tb(i)) * real(Y(tb(i), fb(i))) + abs(V(fb(i)) * Y(tb(i), fb(i)))*cos(angle(Y(tb(i), fb(i))) + delta(fb(i)) - delta(tb(i)));
end

%% dQij/ddelta
for i = 1:nbranch
    if fb(i)~=1
        H11(i,fb(i)-1) = abs(V(fb(i)) * V(tb(i)) * Y(fb(i), tb(i))) * cos(angle(Y(fb(i), tb(i))) + delta(tb(i)) - delta(fb(i)));
    end
    if   tb(i)~=1
        H11(i,tb(i)-1) = -abs(V(fb(i)) * V(tb(i)) * Y(fb(i), tb(i))) * cos(angle(Y(fb(i), tb(i))) + delta(tb(i)) - delta(fb(i)));
    end
end

%% dQij/dV
for i = 1:nbranch
    H12(i,tb(i)) = -abs(V(fb(i)) * Y(fb(i), tb(i))) * sin(angle(Y(fb(i), tb(i))) + delta(tb(i)) - delta(fb(i)));
    H12(i,fb(i)) = -2 * V(fb(i)) * (0.5*abs(b(i)) - imag(Y(fb(i), tb(i)))) - abs(V(tb(i)) * Y(fb(i), tb(i))) * sin(angle(Y(fb(i), tb(i))) + delta(tb(i)) - delta(fb(i)));
end

%% dQji/ddelta
for i = 1:nbranch
    if tb(i)~=1
        H13(i,tb(i)-1) = abs(V(tb(i)) * V(fb(i)) * Y(tb(i), fb(i))) * cos(angle(Y(tb(i), fb(i))) + delta(fb(i)) - delta(tb(i)));
    end
    if   fb(i)~=1
        H13(i,fb(i)-1) = -abs(V(tb(i)) * V(fb(i)) * Y(tb(i), fb(i))) * cos(angle(Y(tb(i), fb(i))) + delta(fb(i)) - delta(tb(i)));
    end
end

%% dQji/dV
for i = 1:nbranch
    H14(i,fb(i)) = -abs(V(tb(i)) * Y(tb(i),fb(i))) * sin(angle(Y(tb(i), fb(i))) + delta(fb(i)) - delta(tb(i)));
    H14(i,tb(i)) = -2 * V(tb(i)) * (0.5 * abs(b(i)) - imag(Y(tb(i), fb(i)))) - abs(V(fb(i)) * Y(tb(i), fb(i))) * sin(angle(Y(tb(i), fb(i))) + delta(fb(i)) - delta(tb(i)));
end

H=[ H1, H2; ...
    H3, H4; ...
    H5, H6; ...
    H7, H8; ...
    H9, H10; ...
    H11, H12; ...
    H13, H14];

end