

function Dist = src_calculaMatrizDistancias(P)

[npuntos nada] = size(P);
Dist = zeros(npuntos, npuntos);

for x=1:npuntos
    for y=1:npuntos
        resta = P(x,:) - P(y,:);
        mult = resta .* resta;
        Dist(y,x) = sqrt( mult(1) + mult(2) );
    end
end


