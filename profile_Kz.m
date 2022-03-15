function [Kz, Kz0, zh] = profile_Kz(z, Kz0_temp, n, zm, zh_temp)

    if 2*(1+n)*zm/(2+n)<zh_temp
        display('Error: check zh and zm values')
    end

    Kz = zeros(length(z),1); 
    
    if n<=1e-10 % n==0
        zh = zm;
        % Parámetro a queda libre
        % Region I
        a = -1/zm^2;
        KzI = Kz0_temp*(1 - 2*a*zm*z + a*z.^2);
        Kz(z>zh) = KzI(z>zh);

        % Region II
        Kz0 = 2*Kz0_temp;   
        Kz(z<=zh) = Kz0;

    else
        Kz0 = Kz0_temp;
        zh = zh_temp;

        % Region I        
        a = 1/(2*zh/n*(zm-zh)-zh*(zh-2*zm));
        KzI = Kz0*(1 - 2*a*zm*z + a*z.^2);
        Kz(z>zh) = KzI(z>zh);

        % Region II
        e = a*2*(zm-zh)*zh/n;
        KzII = Kz0*e*abs(z/zh).^(-n);   
        Kz(z<=zh) = KzII(z<=zh);
    end

end
