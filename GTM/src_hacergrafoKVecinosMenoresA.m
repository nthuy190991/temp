% P: Matriz de r x 2 con las coordenadas (x,y) de los puntos
% conectividad: entero que indica el grado que tendra cada vertice
% Dist: Matriz cuadrada de distancias entre todos los puntos
function Adj = src_hacergrafoKVecinosMenoresA(P,conectividad, Dist, mediana)

[npuntos nada] = size(P);

Adj  = zeros(npuntos, npuntos);
P = [];

for col=1:npuntos
    [VAL PERM] = sort(Dist(:,col),1,'ascend'); %VAL es el vector ordenado, PERM son los indices ordenados
    punto = 2;
    num_conecciones = 0; 
    %while num_conecciones <= conectividad && punto <= npuntos
    %P = [];
    while num_conecciones < conectividad && punto <= npuntos
        %if col == PERM(punto)
        %if VAL(punto) == 0
            %punto = punto + 1;
        if (VAL(punto) <= mediana)
            Adj(col,PERM(punto)) = 1;
                Adj(PERM(punto),col) = 1;
            num_conecciones = num_conecciones + 1;
            punto = punto + 1;
        else
            punto = npuntos+1; %si es mayor que la mediana, ya no sigue buscando
        end 
    end
            
    %if num_conecciones<=conectividad  %Desconecta a todo el vertice
    if num_conecciones<conectividad  %Desconecta a todo el vertice
        P = [P col]; %guarda el indice del vertice a desconectar
    end
end

% Desconecta a todos los vertices que no cumplieron
[r c] = size(P);

for j=1:c
    for k=1:npuntos
       Adj(k, P(j)) = 0;
       Adj(P(j), k) = 0;
    end
end

%Limpia la diagonal
for j=1:npuntos
    Adj(j,j) = 0;
end
