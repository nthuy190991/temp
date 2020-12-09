
function [NP1, NP2] = src_EliminaPuntosRepetidos(P1,P2) 

% Elimina las relaciones 1 a muchos
[nv nc] = size(P1);
NP1=[];
NP2=[];
for i=1:nv
    igual = 0;
    for j=i+1:nv  
        if (P1(i,1)==P1(j,1) && P1(i,2)==P1(j,2)) || (P2(i,1)==P2(j,1) && P2(i,2)==P2(j,2))
            igual = 1;
        end
    end
    if igual==0
        NP1 = [NP1; P1(i,1:2)];
        NP2 = [NP2; P2(i,1:2)];
    end
end
