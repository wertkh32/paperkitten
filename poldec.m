function [R S] = poldec(F)
[A B C] = svd(F);

%D = [1 0 0;0 1 0;0 0 det(C * A')];
if det(F) < 0
    [~, h] = min(diag(B));
    C(:,h) = -C(:,h);
end

R = A * C';
S = C * B * C';
S = (S + S')/2;
end
