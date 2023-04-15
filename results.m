%% "UKF Parameter Tuning for Local Variation Smoothing"
%% K.Nielsen, 2021

x = linspace(-4,4,1000)';
a = 4*pi;
y = 0.1*x.^3 + 0.4*sin(a*x);

% x ~ N(2.2, 0.25) ratio W/sigma faible
xp = 2.2 + sqrt(0.25)*randn(100,1);
mux = mean(xp);
Px = mean((xp-mean(xp)*ones(100,1)).^2);

yp = 0.1*xp.^3 + 0.4*sin(a*xp);
muy = mean(yp);
yp = muy + sqrt(0.01)*randn(100,1);

% x ~ N(-2.2, 0.0025) ratio W/sigma élevée
XP = -2.2 + sqrt(0.0025)*randn(100,1);
MUX = mean(XP);
PX = mean((XP-mean(XP)*ones(100,1)).^2);

YP = 0.1*XP.^3 + 0.4*sin(a*XP);
MUY = mean(YP);
YP = MUY + sqrt(0.01)*randn(100,1);

% représentation graphique de la fonction
figure(1);
plot(x,y,'k-');
xlabel("x");
ylabel("y");
hold on;
plot(xp,yp,'bx');
plot(mux,muy,'ro','LineWidth',1.5);
plot(XP,YP,'bx');
plot(MUX,MUY,'ro','LineWidth',1.5);

%-----------------------------------------------------------
alpha = 1;
beta = 0;
kappa = 0;

[muz,Pz] = sut(2.2,0.25,alpha,beta,kappa,a);

%plot(2.2,muz,"cs");

[MC_muz,MC_Pz] = MC(2.2,0.25,alpha,beta,kappa,a,100000);

plot(2.2,MC_muz,"cs");

muz
Pz
disp("-----")
MC_muz
MC_Pz

%%
%---- Figure 2 ------
[tab_muz,tab_Pz,param] = data_fig2(2.2,0.25,4*pi);
figure(2); hold on;
plot(param,tab_muz,'b-', ...
     param,tab_Pz,'r-', ...
     param,MC_muz*ones(size(param)),'b--', ...
     param,MC_Pz*ones(size(param)),'r--');
xlabel('\alpha√(1+\kappa)');
legend('\mu_z','P_z','\mu_z (MC)', 'P_z (MC)','Location','SouthEast');
%%
%----- Figure 3 ----
[tab_Pz_beta,param] = data_fig3(2.2,0.25,4*pi);
figure(3);
semilogy(param,tab_Pz_beta);
hold on;
plot(param,MC_Pz*ones(size(param)),'b--');
xlabel('\alpha√(1+\kappa)');
ylabel('P_z');
legend('\beta = 0','\beta = 0.01','\beta = 0.027','\beta = 0.072','\beta = 0.19','\beta = 0.52','\beta = 1.39','\beta = 3.73','\beta = 10');
%%
%---- Figure 4 -----
[tab_muz_ratio, tab_Pz_ratio, param] = data_fig4(2.2,4*pi);
figure(4);
plot(param,tab_muz_ratio);
hold on;
plot(param,MC_muz*ones(size(param)),'b--');
xlabel('\alpha√(1+\kappa)');
ylabel('\mu_z');
legend('W/\sigma = 50','W/\sigma = 5','W/\sigma = 2.5','W/\sigma = 1.67','W/\sigma = 1.25','MC','Location','SouthEast');
%% Figure 5
N = 100;
x0 = 2.2;
P0 = 0.25;
x = x0 + sqrt(P0)*randn(N,1);
xinit = x-x0*ones(N,1);

UT1 = zeros(N,1);
UT2 = zeros(N,1);
CT = zeros(N,1);
PSE = [UT1,UT2,CT];
ALPHA = [1;1e-3;1];
BETA = [0;2;0];
KAPPA = [2;0;0];
for k = 1:2:3
    for i = 1:N
        x = x0 + sqrt(P0)*randn(N,1);
        [PS,unused] = ukf(mean(x),mean((x-mean(x)*ones(N,1)).^2),ALPHA(k),BETA(k),KAPPA(k),4*pi);
        [MC_PS,unused] = MC(x0,P0,ALPHA(k),BETA(k),KAPPA(k),4*pi,N);
        PSE(i,k) = PS - MC_PS;
    end
end

g1 = repmat({'xinit'},size(xinit));
g2 = repmat({'UT1'},size(xinit));
g3 = repmat({'UT2'},size(xinit));
g4 = repmat({'CT'},size(xinit));

X = [xinit;PSE(:,1);PSE(:,2);PSE(:,3)];
G = [g1;g2;g3;g4];

boxplot(X,G);

%% _____________________________________________________________
function [tab_muz,tab_Pz,param] = data_fig4(x,a)
    RATIO = [50;5;2.5;1.67;1.25];
    SIGMA = 0.5./RATIO;
    tab_muz = zeros(200,5);
    tab_Pz = zeros(200,5);
    param = zeros(200,1);
    ALPHA = linspace(0.01,1,200);
    beta = 0;
    kappa = 0;
    for i = 1:5
        P = SIGMA(i)^2;
        for j = 1:length(ALPHA)
            param(j) = ALPHA(j)*sqrt(1+kappa);
            [tab_muz(j,i),tab_Pz(j,i)] = sut(x,P,ALPHA(j),beta,kappa,a);
        end
    end
end

function [tab_Pz_beta,param] = data_fig3(x,P,a)
    tab_Pz_beta = zeros(200,9);
    BETA = [0,0.01,0.027,0.072,0.19,0.52,1.39,3.73,10];

    for i = 1:9
        tab1_muz = zeros(100,1);
        tab1_Pz = zeros(100,1);
        param1 = zeros(100,1);
        ALPHA = linspace(0.01,1,100);
        for j = 1:length(ALPHA)
            kappa = 0;
            param1(j) = ALPHA(j)*sqrt(1+kappa);
            [tab1_muz(j),tab1_Pz(j)] = sut(x,P,ALPHA(j),BETA(i),kappa,a);
        end
        tab2_muz = zeros(100,1);
        tab2_Pz = zeros(100,1);
        param2 = zeros(100,1);
        KAPPA = linspace(0,10,100);
        for j = 1:length(KAPPA)
            alpha = 1;
            param2(j) = alpha*sqrt(1+KAPPA(j));
            [tab2_muz(j),tab2_Pz(j)] = sut(x,P,alpha,BETA(i),KAPPA(j),a);
        end
        tab_muz = [tab1_muz;tab2_muz];
        tab_Pz = [tab1_Pz;tab2_Pz];
        param = [param1;param2]; 

        tab_Pz_beta(:,i) = tab_Pz;
    end
end


function [tab_muz,tab_Pz,param] = data_fig2(x,P,a)
    tab1_muz = zeros(100,1);
    tab1_Pz = zeros(100,1);
    param1 = zeros(100,1);    
    ALPHA = linspace(0.01,1,100);
    for i = 1:length(ALPHA)
        beta = 0;
        kappa = 0;
        param1(i) = ALPHA(i)*sqrt(1+kappa);
        [tab1_muz(i),tab1_Pz(i)] = sut(x,P,ALPHA(i),beta,kappa,a);
    end

    tab2_muz = zeros(100,1);
    tab2_Pz = zeros(100,1);
    param2 = zeros(100,1);
    KAPPA = linspace(0,10,100);
    for i = 1:length(KAPPA)
        alpha = 1;
        beta = 0;
        param2(i) = alpha*sqrt(1+KAPPA(i));
        [tab2_muz(i),tab2_Pz(i)] = sut(x,P,alpha,beta,KAPPA(i),a);
    end

    tab_muz = [tab1_muz;tab2_muz];
    tab_Pz = [tab1_Pz;tab2_Pz];
    param = [param1;param2];
end

%-----------------

function [MC_muz,MC_Pz] = MC(meanx,covx,alpha,beta,kappa,a,N)
    tab_muz = zeros(N,1);
    tab_Pz = zeros(N,1);
    for i = 1:N
        xp = meanx + sqrt(covx)*randn(100,1);
        mux = mean(xp);
        Px = mean((xp-mean(mux)*ones(100,1)).^2);
        [tab_muz(i),tab_Pz(i)] = sut(mux,Px,alpha,beta,kappa,a);
    end
    MC_muz = mean(tab_muz);
    MC_Pz = mean(tab_Pz);
end

%-----------------

function [y,Py] = sut(x,P,alpha,beta,kappa,a)
    nx = length(x);
    xsp = zeros(2*nx+1,nx);
    Wsp = zeros(2*nx+1,nx);
    Wc = zeros(2*nx+1,nx);

    lambda = alpha^2*(nx+kappa)-nx;

    % calcul des poids
    Wsp(1) = lambda/(nx+lambda);
    Wc(1) = lambda/(nx+lambda) + (1-alpha^2 + beta);
    for i = 2:2*nx+1
        Wsp(i) = 1/(2*(nx+lambda));
        Wc(i) = 1/(2*(nx+lambda));
    end

    % calcul des SP
    delta = chol((nx+lambda)*P)';
    xsp(1) = x;
    for i = 1:nx
        xsp(i+1) = x + delta;
        xsp(i+nx+1) = x - delta;
    end

    % propagation des SP
    ysp = zeros(2*nx+1,nx);
    for i = 1:2*nx+1
        ysp(i) = 0.1*xsp(i)^3 + 0.4*sin(a*xsp(i));
    end
    %plot(xsp,ysp,'mo','LineWidth',2)

    % moyenne et covariance des SP propagés
    y = 0;
    for i = 1:2*nx+1
        y = y + Wsp(i)*ysp(i);
    end
    Py = 0;
    for i = 1:2*nx+1
        Py = Py + Wc(i)*(ysp(i)-y)*(ysp(i)-y)';
    end
    %Py = Py + (1-alpha^2 + beta)*(ysp(1)-y)*(ysp(1)-y)';
end


function [post_x, post_P] = ukf(x,P,alpha,beta,kappa,a)
    nx = length(x);
    xsp = zeros(2*nx+1,nx);
    Wsp = zeros(2*nx+1,nx);
    Wc = zeros(2*nx+1,nx);

    lambda = alpha^2*(nx+kappa)-nx;

    % calcul des poids
    Wsp(1) = lambda/(nx+lambda);
    Wc(1) = lambda/(nx+lambda) + (1-alpha^2 + beta);
    for i = 2:2*nx+1
        Wsp(i) = 1/(2*(nx+lambda));
        Wc(i) = 1/(2*(nx+lambda));
    end

    % calcul des SP
    delta = chol((nx+lambda)*P)';
    xsp(1) = x;
    for i = 1:nx
        xsp(i+1) = x + delta;
        xsp(i+nx+1) = x - delta;
    end

    % propagation des SP dans la fonction de prédiction
    pred_xsp = zeros(2*nx+1,nx);
    for i = 1:2*nx+1
        pred_xsp(i) = 0.1*xsp(i)^3 + 0.4*sin(a*xsp(i));
    end
    %plot(xsp,ysp,'mo','LineWidth',2)

    % moyenne et covariance des SP prédits
    pred_x = 0;
    for i = 1:2*nx+1
        pred_x = pred_x + Wsp(i)*pred_xsp(i);
    end 
    pred_P = 0;
    for i = 1:2*nx+1
        pred_P = pred_P + Wc(i)*(pred_xsp(i)-pred_x)*(pred_xsp(i)-pred_x)';
    end

    %-- passage des SP prédits dans la fonction d'observation
    ysp = zeros(2*nx+1,nx);
    for i = 1:2*nx+1
        ysp(i) = 0.1*pred_xsp(i)^3 + 0.4*sin(a*pred_xsp(i));
    end

    % moyenne et covariance des observations prédites
    pred_y = 0;
    for i = 1:2*nx+1
        pred_y = pred_y + Wsp(i)*ysp(i);
    end 
    Py = 0;
    for i = 1:2*nx+1
        Py = Py + Wc(i)*(ysp(i)-pred_y)*(ysp(i)-pred_y)';
    end

    % cross-covariance
    Pxy = 0; % zeros(nx,nx) pour la prochaine fois
    for i = 1:2*nx+1
        Pxy = Pxy + Wc(i)*(pred_xsp(i)-pred_x)*(ysp(i)-pred_y)';
    end

    % Gain de Kalman
    K = Pxy/Py;

    % mesure y
    ymes = 0.1*x^3 + 0.4*sin(a*x) + sqrt(0.01)*randn(1,1);

    % correction
    post_x = pred_x + K*(pred_y - ymes);
    post_P = pred_P - K*Py*K';
end


function val = h1(muxhat,Pxhat,muxgt)
    nx = length(muxhat);
    val = 1/2*(nx*log(2*pi)+log(det(Pxhat)) + (muxhat-muxgt)*inv(Pxhat)*(muxhat-muxgt)');
end

function val = h2(muyhat,Pyhat,muygt)
    ny = length(muyhat);
    val = 1/2*(ny*log(2*pi)+log(det(Pyhat)) + (muyhat-muygt)*inv(Pyhat)*(muyhat-muygt)');
end

function val = costfun(theta)
    x = 2.2 + sqrt(0.25)*randn(1,1);
    muxgt = 0.1*x^3 + 0.4*sin(4*pi*x);
    [muxhat,Pxhat] = ukf(mux0,Px0,theta(1),theta(2),theta(3),4*pi);
    score = h1(muxhat,Pxhat,muxgt);
    for j = 2:100
        x = 2.2 + sqrt(0.25)*randn(1,1);
        muxgt = 0.1*x^3 + 0.4*sin(4*pi*x);
        [muxhat,Pxhat] = ukf(mux0,Px0,theta(1),theta(2),theta(3),4*pi);
        score = score + h1(muxhat,Pxhat,muxgt);
    end

end