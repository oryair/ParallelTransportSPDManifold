function dist = RiemannianDist(mP1, mP2, p)

   if nargin < 3
       p = 2;
   end

%     mP   = mP1^-1 * mP2;
%     vLam = eig(mP);
    
    vLam = eig(mP2, mP1);
    
    if p == 1 %-- just for speed
        dist = sum( abs(log(vLam)) );
    else
        dist = sum( abs(log(vLam)).^p ).^(1/p);
    end
    
end