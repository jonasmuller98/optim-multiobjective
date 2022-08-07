%% Bee Colony Optimization
% Jonas Muller Gonçalves
clear all;clc;
% Variáveis Auxiliares
parada=0; iter=0; iparada=0; itparada=0;

%% Dados da Função de Otimização
% f = @(x) 100*(x(2)-x(1).^2).^2+(1-x(1)).^2; % Função Objetivo (**modificar o plot no final do programa**)
dominio=[-2 4;-2 4];                     % Domínino da primeira variável de projeto
% dominio(2,:) = [-2 4;-2 4];  % Domínino da segunda variável de projeto
tol = 1e-6;     % Tolerância da diferença da f(n-1) e f(n)
n_abelhas = 5;   % Nº de operárias na colméia
n_variaveis = 2;  % Variáveis de projeto
n_seguidoras = 2; % Número de seguidoras da colméia
iter_final = 100;    % Número de iterações com f(n-1)-f(n) < tol
ipart=n_abelhas;    % auxiliar.

for p1=0:0.1:1
    f1 = @(x) x(1).^2;
    f2 = @(x) (x(1)-2).^2;
    
    f = @(x) p1*(x(1).^2) + (1-p1)*(x(1)-2).^2;
    
    %% Inicialização das Operárias Iniciais
    for icol = 1:n_variaveis
        for ivar = 1:n_abelhas
            operarias_novas(ivar,icol) = dominio(icol,1) + rand*(dominio(icol,1)-dominio(icol,2));
        end
    end
    
    %% Loop Principal
    while iparada==0
        iter=iter+1;
        operarias_inicial=operarias_novas;
        
        % Confere se as operárias estão no domínio
        for ipart=1:n_abelhas
            for ivariab=1:n_variaveis
                if operarias_novas(ipart,ivariab) < dominio(ivariab,1)
                    operarias_inicial(ipart,ivariab) = rand*dominio(ivariab,1);  % Induz o retorno da partícula que passou
                elseif operarias_novas(ipart,ivariab) >  dominio(ivariab,2)
                    operarias_inicial(ipart,ivariab) = rand*dominio(ivariab,2);  % Induz o retorno da partícula que passou
                end
            end
        end
        
        %% Avalia a função nas Operárias Iniciais
        for ifunc=1:n_abelhas
            resp_oper_inicial(ifunc,1) = f(operarias_inicial(ifunc,:));
        end
        
        %% Lança as Seguidoras Iniciais
        for i=1:n_abelhas
            for n=1:n_seguidoras
                for j=1:n_variaveis
                    seguidoras_inicial(n,j,i) = operarias_inicial(i,j)+rand;
                end
            end
        end
        
        % Confere se as seguidores estão no domínio
        for ipart=1:n_abelhas
            for isegd=1:n_seguidoras
                for ivariab=1:n_variaveis
                    if seguidoras_inicial(isegd,ivariab,ipart) < dominio(ivariab,1)
                        seguidoras_inicial(isegd,ivariab,ipart) = rand*dominio(ivariab,1);  % Induz o retorno da partícula que passou
                    elseif seguidoras_inicial(isegd,ivariab,ipart) >  dominio(ivariab,2)
                        seguidoras_inicial(isegd,ivariab,ipart) = rand*dominio(ivariab,2);  % Induz o retorno da partícula que passou
                    end
                end
            end
        end
        
        % Avalia a função nas Seguidoras Iniciais
        for iabelha=1:n_abelhas
            for ifunc = 1:n_seguidoras
                resp_seg_inicial(ifunc,iabelha) = f(seguidoras_inicial(ifunc,:,iabelha));
            end
        end
        
        %% Comparação das funções entre seguidoras e operárias
        for i=1:n_abelhas
            for j=1:n_seguidoras
                if resp_oper_inicial(i,1) < min(resp_seg_inicial(:,i))
                    respostas(i,1) = resp_oper_inicial(i,1);
                    respostas(i,2) = i;
                    respostas(i,3) = 1;
                    contador(i,1) = 1; % contador verifica se foi ativado resp oper
                    
                elseif resp_oper_inicial(i,1) > min(resp_seg_inicial(:,i))
                    respostas(i,1) = min(resp_seg_inicial(:,i));
                    [~,respostas(i,2)] = min(resp_seg_inicial(:,i));
                    respostas(i,3) = i;
                    contador(i,1) = 0; % contador verifica se foi ativado resp seg
                end
            end
        end
        
        %% Cálculo das Probabilidades
        for iab=1:n_abelhas
            Pi_aux(iab,1) =  1-(respostas(iab,1)/sum(respostas(:,1))); % probabilidade auxiliar
        end
        
        Pi_Roleta(1,1) = 0; % seta a primeira probabilidade como 0
        
        for i=1:n_abelhas
            Pi_Roleta(i+1,1) = sum(Pi_aux(1:i,1));  % calcula as próximas probabilidades de 0 a 1
        end
        
        %     Gira a Roleta e atribui o Valor das Posições de Acordo com as Probabilidades
        
        for j=1:n_abelhas
            aleatorio = rand;
            for i=1:n_abelhas
                if Pi_Roleta(i,1) < aleatorio && aleatorio <= Pi_Roleta(i+1,1)
                    a=1; armazena(j,1) = Pi_Roleta(i,1);
                    armazenalinha1(j,1) = i;
                    contador2(j,1) = contador(i,1);
                    colunas(j,1) = respostas(armazenalinha1(j,1),3);
                    linhas(j,1) = respostas(armazenalinha1(j,1),2);
                    
                    if contador2(j,1) == 1
                        operarias_novas(j,:) = operarias_inicial(linhas(j,1),:);
                    elseif contador2(j,1) == 0
                        operarias_novas(j,:) = seguidoras_inicial(linhas(j,1),:,colunas(j,1));
                    end
                end
            end
        end
        
        for iconfer = 1:n_abelhas
            armazenaf_aux(iconfer,1) = f(operarias_novas(iconfer,:));
            [ordenados,indice] = sort(armazenaf_aux);
            armazenaf(iter,1) = min(armazenaf_aux);
        end
        operaria = operarias_novas(indice(1,1),:);
        
        %% Critério de Parada
        if iter ~= 1
            if armazenaf(iter,1) - armazenaf(iter-1,1) < tol
                itparada = itparada + 1;
                if itparada == iter_final
                    iparada=1;
                end
            else
                itparada=0;
            end
        end
        
        %% Lançando exploradoras em iterações aleatórias
        for icol = 1:n_variaveis
            if iter == randi([iter iter+10])
                operarias_novas(indice(end,1),icol) = dominio(icol,1) + rand*(dominio(icol,1)-dominio(icol,2));
            end
        end
        operaria2(iter,:) = operaria;
        for ib=1:iter
            fs1(ib,1) = f1(operaria2(ib,:));
        end
        for ic=1:iter
            fs2(ic,1) = f2(operaria2(ic,:));
        end
    end
    
    
    
    scatter(min(fs1),min(fs2));
    hold on;
    
    
end