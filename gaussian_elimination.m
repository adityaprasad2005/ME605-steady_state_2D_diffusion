function x = gaussian_elimination(A, b)
    % Gaussian elimination with pivoting

    [n, ~] = size(A);  % Number of equations

    % Forward elimination
    for i = 1:1:n-1
        for j = i+1:1:n
        factor = A(j,i)/A(i,i);
            for k = i:1:n
                A(j,k) = A(j,k) - factor*(A(i,k));
            end
            b(j) = b(j) - factor*(b(i));    
        end
    end

    % Back substitution:
    x = zeros(n,1);

    for i = n:-1:1
        s = 0;
        for j = i+1:1:n
            s = s + A(i,j)*x(j);
        end
        x(i) = (b(i) - s)/A(i,i);
    end
end
