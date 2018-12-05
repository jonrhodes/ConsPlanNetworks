function BR=buildSpecies(n,s,p,q)
%BR = a site x species matrix of the community
%n = the number of sites
%s = maximum number of species
%p = expected probability that a species occurs at a site
%q = nestedness with values between 0 and 1
%(0 = maximum nestedness to 0.5 = random to 1.0 = maximum endemism)
%note that in the paper we call 1-q nestedness so that maximum nestedness =
%1 and minimum nestedness = 0 (although can't equal zero - see note
%below)
    SP=binornd(s,p,n,1);
    [SP IX]=sort(SP);
    BR=zeros(n,s);
    SubBR=BR(n,:);
    SubBR(randsample(s,SP(n)))=1;
    BR(n,:)=SubBR;
    SEL=SubBR;
    w=SEL-SEL*q+abs(SEL-1)*q;
    for i=1:n-1
        SubBR=BR(n-i,:);
        w_temp=w;
        if SP(n-i) > 0
            for j=1:SP(n-i)
                if all(w_temp == 0)
                    %this causes issues when q = 1 as it makes the distribution more random like, so q can approach 1, but not = 1
                    Num=randsample(s,1,true,abs(SubBR-1));
                else
                   Num=randsample(s,1,true,w_temp);
                end
                w_temp(Num)=0;
                SubBR(Num)=1;
            end
            BR(n-i,:)=SubBR;
            SEL=SubBR;
            w=w.*(SEL-SEL*q+abs(SEL-1)*q);
        end
    end
    BR=BR(IX,:);
end
