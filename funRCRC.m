function zdata = funRCRC(p,w)
    % this function returns a nx2 matrix response of an equiavlent EIS circuit
    % 
    % Inputs:
    % p = parameters to fit. [r-sol R1 C1 R2 C2]
    % w = measurement frequencies (rad/s)
    %
    % Output:
    % zdata = nx2 measured impedance at each frequency. Data is in the form
    % x+i*y. where the data in the first column is the "real" term x, and
    % the data in the second column is the "imaginary" term y.

    % convert the input signal to the correct "direction"
    [n,m] = size(w);
    if n==1 && m~=1
        w = w';
    end
    
%     global avg_alpha
%     alpha = avg_alpha;
    alpha = 1;
    %p = exp(p);  % used to normalize the magnitude of the contribution between the capacitance and resistance terms
    p = 10.^(p);  % used to normalize the magnitude of the contribution between the capacitance and resistance terms
    Z = p(1) + p(2)./(1+(p(2).*p(3).*1i.*w).^alpha) + p(4)./(1+(p(4).*p(5).*1i.*w).^alpha);
    
    zdata(:,1) = real(Z);
    zdata(:,2) = imag(Z);

end