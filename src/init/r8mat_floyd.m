function B = r8mat_floyd(A)
%R8MAT_FLOYD Shortest-path distances for a dense adjacency matrix.
% Missing edges should be Inf. This local implementation is included so the
% public paper repo has no third-party shortest-path dependency.
    B = A;
    n = size(A, 1);
    for k = 1:n
        B = min(B, B(:, k) + B(k, :));
    end
end
