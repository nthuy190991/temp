function [P_final, P2_final, MKNNG_0, MKNNG2_0, MKNNG, MKNNG2] = GTM_algo_original(P_init, P2_init, K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source code of Wendy Aguilar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% La diferencia con la propuesta 1 es que ahora crea el grafo conectando
% un vertice con sus k vecinos mas cercanos siempre y cuando la distancia
% sea menor a la mediana. Si no encuentra k vecinos que cumplan con esta
% restriccion lo toma como outlier.

%% Input
conectividad = K;
P1 = fliplr(P_init);
P2 = fliplr(P2_init);

    d1 = distmat(P_init);
    d2 = distmat(P2_init);

    median1 = median(d1(:));
    median2 = median(d2(:));

    [MKNNG_0]  = MedianKNearestNeighborGraph(d1,K,median1);
    [MKNNG2_0] = MedianKNearestNeighborGraph(d2,K,median2);

%% Body program
[NP1, NP2] = src_EliminaPuntosRepetidos(P1,P2);

[nv nada] = size(NP1);
match = eye(nv,nv);
% paintMatchVerticalColor(I1,I2,NP1,NP2,match, K);

umbral = 0;
num_iterac = 0;
salida = 0;

t_ini = clock;

Dist1 = src_calculaMatrizDistancias(NP1);
Dist2 = src_calculaMatrizDistancias(NP2);

mediana1 = src_obtenMediana(Dist1);
mediana2 = src_obtenMediana(Dist2);

% size(NP1)
while salida==0
    num_iterac = num_iterac + 1;
    % Crea los dos grafos 
    G1 = src_hacergrafoKVecinosMenoresA(NP1, conectividad, Dist1, mediana1);
    G2 = src_hacergrafoKVecinosMenoresA(NP2, conectividad, Dist2, mediana2);
 
    [m m] = size(G1);

    % Crea la matriz de las estructuras en comun
    comun = abs(G1-G2);
        
    %quita el vertice con mayor cant de 1's en su columna, es outlier
    suma = sum(comun);
    if sum(suma) > umbral
        indices = find(suma==max(suma));
        indi = indices(1)-1;
        indf = indices(1)+1;
        NP1 = [NP1(1:indi,1:2) ; NP1(indf:m,1:2)];
        NP2 = [NP2(1:indi,1:2) ; NP2(indf:m,1:2)];
        Dist1 = [Dist1(1:indi,1:indi) Dist1(1:indi,indf:m); Dist1(indf:m,1:indi) Dist1(indf:m,indf:m)];
        Dist2 = [Dist2(1:indi,1:indi) Dist2(1:indi,indf:m); Dist2(indf:m,1:indi) Dist2(indf:m,indf:m)];
    else
        salida = 1; %termina
    end
end

t_fin = clock;
tiempo = etime(t_fin,t_ini);


%% Output
    P_final  = fliplr(NP1);
    P2_final = fliplr(NP2);

    d1 = distmat(P_final);
    d2 = distmat(P2_final);

    median1 = median(d1(:));
    median2 = median(d2(:));

    [MKNNG]  = MedianKNearestNeighborGraph(d1,K,median1);
    [MKNNG2] = MedianKNearestNeighborGraph(d2,K,median2);
