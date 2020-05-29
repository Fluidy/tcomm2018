function [ q ] = Syn_preference( p, s, theta )
Nf = length(p);
Nu = length(s);
q = zeros(Nu, Nf);
pp = p;
pp0 = pp;

availabe_user = 1:Nu;
for i_ue = 1 : Nu
    ii = randperm(length(availabe_user), 1);
    ue = availabe_user(ii);
    availabe_user(ii) = [];
  
    rho = min(pp/s(ue), 1);
    g = zeros(1, Nf);
    aa = 1;
    availabe_file = 1:Nf;

    while aa > 1e-8 && ~isempty(availabe_file)
        ff = randperm(length(availabe_file), 1);
        f = availabe_file(ff);
        
        scale = (aa/sum(pp0(availabe_file))).^theta; 
        
        x = theta*pp0(f);
        y = theta*pp0(f) + (1-theta)*rho(f);       
        temp = min(min(aa, scale*(x + (y-x)*rand)), rho(f));
        
        g(f) = temp;
        
        aa = aa - g(f);
        availabe_file(ff) = [];

    end

    if aa > 1e-8 
        not_full_file = find(rho - g > 0);
        while ~isempty(not_full_file) && aa > 1e-8
            ff = randperm(length(not_full_file), 1);
            f = not_full_file(ff);
            temp = min(g(f) + aa, rho(f));
            aa = aa - (temp - g(f));
            g(f) = temp;
            if g(f) == rho(f)
                not_full_file(ff) = [];
            end
        end
    end

    q(ue,:) = g;
    pp = max( pp - s(ue).*q(ue,:), 0);
    
    if min(pp) < 0
        warning('!!');
    end
    
    pp0 = pp/sum(pp);
end

for ue = 1 : Nu
    q(ue, :) = q(ue, :)*s(ue);
end
