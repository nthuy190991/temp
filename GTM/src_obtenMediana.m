
function Mediana = src_obtenMediana(Distancias)

%Coloca todos los elementos en un vector

[npuntos nada] = size(Distancias);
M = [];

for i=1:npuntos
   M = [M Distancias(i,:)]; 
end

Mediana = median(M);





