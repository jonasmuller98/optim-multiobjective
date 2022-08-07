%% Particle Swarm Optimization
% Jonas Müller Gonçalves
clear all;clc;
% Não tenho absoluta certeza a respeito da restrição. A função implementada
% aqui é Mishra's Bird function - Constrained presente no Wikipedia

%% Inicialização das Variáveis
population = 50; npart=population;   % População total de Partículas
variables = 2; nvar=variables;     % Variáveis de Projeto
lambda1 = 2;    % Lambda 1
lambda2 = 2;    % Lambda 2
omega = 0.5;    % Propriedade de Inércia
dominio=[-2 4;-2 4]; % Domínio das nvar variaveis
tol = 1e-8;     % Diferença entre o xgbest anterior e o novo xgbest - erro
erro_contador = 0;  % contador de iteracoes com o verificador de erro ativado
nerros = 200;   % Número máximo de iterações com o erro ativado
% rng(2);

for p1=0:0.1:1 
f1 = @(x) x(1).^2;
f2 = @(x) (x(1)-2).^2;
f = @(x) p1*x(1).^2 + (1-p1)*(x(1)-2).^2;

%% Inicialização dos Vetores de Projeto
% Distribuindo as Partículas
for icol = 1:nvar
    for ivar = 1:npart
        position(ivar,icol) = dominio(icol,1)+rand*(dominio(icol,1)-dominio(icol,2));
    end
end
ct1=1;ct2=1;    % contadores auxiliares

% Checando se o chute inicial está dentro do domínio
for ipart=1:npart
    for ivariab=1:nvar
        if position(ipart,ivariab) < dominio(ivariab,1)
            fora_dominio_menor(ct1,ivariab) = position(ipart,ivariab); ct1=ct1+1;
            position(ipart,ivariab) = rand*dominio(ivariab,1);  % Induz o retorno da partícula que passou
        elseif position(ipart,ivariab) >  dominio(ivariab,2)
            fora_dominio_maior(ct2,ivariab) = position(ipart,ivariab); ct2=ct2+1;
            position(ipart,ivariab) = rand*dominio(ivariab,2);  % Induz o retorno da partícula que passou
        end
    end
end

% Definindo as Velocidades Iniciais
velocity = zeros(npart,nvar);   % inicializa as velocidades

%% Cálculo Inicial da Melhor Partícula
for ires=1:size(position,1)
    resultados(ires,1) = f(position(ires,:));    % Define os valores de f(x,y)
end

[resg,ord]=sort(resultados);    % Ordena os valores do menor para o maior
for ixgbest=1:nvar
    xgbest(ixgbest) = [position(ord(1,1),ixgbest)];   % Define o melhor global como o primeiro do vetor ordenado
end
melhorglobal=xgbest;    % Cria um valor armazenado do melhor global
ordglobal = ord(1,1);   % Replica qual a linha do melhor global

xlbest = zeros(npart,nvar); % Inicia Vetor de melhor local

% Variáveis Auxiliares
parada=0; iter=0; pos_old=position; v_old=velocity; cont1=1; cont2=1;


%% Otimização
while parada == 0
    iter = iter+1; % Atualiza o Contador de Iterações
    r1=rand; r2=rand; % R's aleatórios
    
    % Cálculo das Novas Velocidades e Posições de cada partícula
    for i=1:npart
        for j=1:nvar
            v(i,j) = omega*v_old(i)+lambda1*r1*(xlbest(i,j) - pos_old(i,j)) + lambda2*r2*(xgbest(j) - pos_old(i,j));
            pos(i,j) = pos_old(i,j) + v(i,j);
        end
    end
    cont1=1;cont2=1;    % Contadores Auxiliares
    
    % Verificação do Domínio das Novas Partículas
    for ipart=1:npart
        for ivariab=1:nvar
            if pos(ipart,ivariab) < dominio(ivariab,1)
                fora_dominio_menor(ct1,ivariab) = pos(ipart,ivariab); ct1=ct1+1;    % armazena os valores que foram menores que o dominio (somente para conferencia)
                pos(ipart,ivariab) = rand*dominio(ivariab,1);  % Induz o retorno da partícula que passou
            elseif pos(ipart,ivariab) >  dominio(ivariab,2)
                fora_dominio_maior(ct2,ivariab) = pos(ipart,ivariab); ct2=ct2+1;    % armazena os valores que foram maiores que o dominio (somente para conferencia)
                pos(ipart,ivariab) = rand*dominio(ivariab,2);  % Induz o retorno da partícula que passou
            end
        end
    end
    
    % Atualiza o valor de xlbest 
    for j=1:nvar
        for i=1:npart
            if f(xlbest(i,:)) > f(pos(i,:))   % Define o melhor global como o primeiro do vetor ordenado
                xlbest(i,:) = pos(i,:);
            else
                
            end
        end
    end
    
    % Armazena os resultados de f(xlbest)
    for ilc=1:npart
        resultados_loc(ilc,1) = f(xlbest(ilc,:));
    end
    
    % Ordena os valores de f(xlbest) do menor para o maior
    [rloc,iloc]=sort(resultados_loc);
    
    % Armazena o melhor x local
    for imelhorl=1:nvar
        melhorlocal(imelhorl) = [xlbest(iloc(1),imelhorl)];
    end
    
    % Confere se o melhor x local é melhor que o melhor x global
    if f(melhorglobal(1,:)) > f(melhorlocal(1,:))
        melhorglobal(1,:) = melhorlocal(1,:);
    end
    
    erro = f(xgbest(1,:)) - f(melhorglobal(1,:));   % Verifica o erro
    xgbest = melhorglobal;    % Atualiza o xgbest 
    v_old = v;  % Atribui a velocidade antiga como a nova velocidade calculada
    pos_old = pos;  % Atribui a posição antiga como a nova posição calculada (já com o domínio verificado)
    melhores(iter,:) = melhorglobal; % Armazena todos melhores globais
    
    % Condição de Parada
    if erro<=tol % se o erro calculado é menor que a tolerancia, atualiza o contador
        erro_contador = erro_contador+1;
        if erro_contador>nerros % se o contador atinge o numero maximo de iteracoes consecutivas com erro<=tol
            parada=1;
        end
    else 
        erro_contador=0;
    end
    

end

%% Pós Processamento

for ires=1:iter
    fmelhores(ires,1) = f(melhores(ires,:));   % Define os valores de f(x,y)
end

for ires=1:iter
    fmelhores1(ires,1) = f1(melhores(ires,:));   % Define os valores de f(x,y)
end

for ires=1:iter
    fmelhores2(ires,1) = f2(melhores(ires,:));   % Define os valores de f(x,y)
end

% figure(1)
% xx=dominio(1,1):.2:dominio(1,2);
% yy=dominio(2,1):.2:dominio(2,2);
% [X,Y]=meshgrid(xx,yy);
% ff=(X - 2).^2;;
% % ff=100*(Y-X.^2).^2+(1-X).^2;
% surf(X,Y,ff); hold on;
% scatter3(melhores(:,1),melhores(:,2),fmelhores(:,1),'r','linewidth',1.4);hold on;
% 
% figure(2);
% scatter(1:iter,fmelhores(:,1));hold on;xlabel('Iteração');ylabel('F Objetivo');grid minor;

figure(3);
scatter(f1(xgbest(1,:)),f2(xgbest(1,:)));hold on;
f(xgbest(1,:))
end

