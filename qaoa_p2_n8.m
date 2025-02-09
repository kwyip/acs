function qaoa_p2_n8(gamma1, gamma2, beta1, beta2)
%function obj = qaoa_p2_n8(gamma1, gamma2, beta1, beta2)
sX = [0 1; 1 0]; sZ = [1 0; 0 -1]; unit = speye(2);
natom = 8; neval = 2^natom; neval2 = neval*2; 

sX_1 = kron(sX, speye(2^(natom-1)));
sX_2 = kron(kron(speye(2),sX),speye(2^(natom-2)));
sX_3 = kron(kron(speye(2^2),sX),speye(2^(natom-3)));
sX_4 = kron(kron(speye(2^3),sX),speye(2^(natom-4)));
sX_5 = kron(kron(speye(2^4),sX),speye(2^(natom-5)));
sX_6 = kron(kron(speye(2^5),sX),speye(2^(natom-6)));
sX_7 = kron(kron(speye(2^6),sX),speye(2^(natom-7)));
sX_8 = kron(speye(2^(natom-1)),sX);

sZsZIIIIII = mkron(sZ,sZ,unit,unit,unit,unit,unit,unit);
IsZsZIIIII = mkron(unit,sZ,sZ,unit,unit,unit,unit,unit);
IIsZsZIIII = mkron(unit,unit,sZ,sZ,unit,unit,unit,unit);
IIIsZsZIII = mkron(unit,unit,unit,sZ,sZ,unit,unit,unit);
IIIIsZsZII = mkron(unit,unit,unit,unit,sZ,sZ,unit,unit);
IIIIIsZsZI = mkron(unit,unit,unit,unit,unit,sZ,sZ,unit);
IIIIIIsZsZ = mkron(unit,unit,unit,unit,unit,unit,sZ,sZ);
sZIIIIIIsZ = mkron(sZ,unit,unit,unit,unit,unit,unit,sZ);

IIIIIIII = speye(2^(natom));

HP = 0.5*(IIIIIIII-sZsZIIIIII) + 0.5*(IIIIIIII-IsZsZIIIII) + 0.5*(IIIIIIII-IIsZsZIIII) + 0.5*(IIIIIIII-IIIsZsZIII) + ...
     0.5*(IIIIIIII-IIIIsZsZII) + 0.5*(IIIIIIII-IIIIIsZsZI) + 0.5*(IIIIIIII-IIIIIIsZsZ) + 0.5*(IIIIIIII-sZIIIIIIsZ);

HB = (sX_1+sX_2+sX_3+sX_4+sX_5+sX_6+sX_7+sX_8);

psi0 = zeros(2^natom,1); 
for i = 1 : 2^natom
    psi0(i) = 1/sqrt(2^natom);
end

% beta1s = -1:0.1:1;
% gamma1s = -1:0.1:1;
% beta2s = -1:0.1:1;
% gamma2s = -1:0.1:1;
beta1s = -1:1.1:1;
gamma1s = -1:1.1:1;
beta2s = -1:1.1:1;
gamma2s = -1:1.1:1;

objs = [];
objs = zeros(length(beta1s),length(gamma1s),length(beta2s),length(gamma2s));

for beta2idx = 1:length(beta2s) 
    for gamma2idx = 1:length(gamma2s)  
        for beta1idx = 1:length(beta1s)
            for gamma1idx = 1:length(gamma1s)
               
                Phi_t = [];
                beta2 = beta2s(beta2idx);
                gamma2 = gamma2s(gamma2idx);
                beta1 = beta1s(beta1idx);
                gamma1 = gamma1s(gamma1idx);
                
                % HP_gamma1 = gamma1.*HP; 
                % HB_beta1 = beta1.*HB;
                % HP_gamma2 = gamma2.*HP; 
                % HB_beta2 = beta2.*HB;
                % 
                % KP_gamma1 = 1*(HP_gamma1^2);
                % KP_gamma1 = [real(KP_gamma1) -imag(KP_gamma1);imag(KP_gamma1) real(KP_gamma1)]; 
                % XP_gamma1 = [zeros(neval2) -eye(neval2);KP_gamma1 zeros(neval2)]; 
                
                XP_gamma1 = toX(gamma1.*HP);

                % KB_beta1 = 1*(HB_beta1^2);
                % KB_beta1 = [real(KB_beta1) -imag(KB_beta1);imag(KB_beta1) real(KB_beta1)]; 
                % XB_beta1 = [zeros(neval2) -eye(neval2);KB_beta1 zeros(neval2)]; 

                XB_beta1 = toX(beta1.*HB);

                % KP_gamma2 = 1*(HP_gamma2^2);
                % KP_gamma2 = [real(KP_gamma2) -imag(KP_gamma2);imag(KP_gamma2) real(KP_gamma2)]; 
                % XP_gamma2 = [zeros(neval2) -eye(neval2);KP_gamma2 zeros(neval2)]; 
                
                XP_gamma2 = toX(gamma2.*HP);

                % KB_beta2 = 1*(HB_beta2^2);
                % KB_beta2 = [real(KB_beta2) -imag(KB_beta2);imag(KB_beta2) real(KB_beta2)]; 
                % XB_beta2 = [zeros(neval2) -eye(neval2);KB_beta2 zeros(neval2)]; 

                XB_beta2 = toX(beta2.*HB);

                % varphi0 = [real(psi0);imag(psi0)];
                % varphidot0_fast = [imag(HP_gamma1)*real(psi0) + real(HP_gamma1)*imag(psi0); -real(HP_gamma1)*real(psi0) + imag(HP_gamma1)*imag(psi0)]; %\varphidot0: initial condition for IVP
                % Phi0 = [varphi0;varphidot0_fast];

                Phi0 = toPhi(gamma1.*HP,psi0);

                if gamma1 ~= 0
                    t = [0 1];
                    [T1,Phi] = ode15s(@(t,Phi)dPhi(t,Phi,XP_gamma1),t,Phi0); %\varphi
                    Phi_t = [Phi_t;Phi];
                    Phi1 = Phi(end, :).';
                else
                    Phi1 = Phi0;
                end
                
                % varphi1 = Phi1(1:neval2);
                % psi1 = varphi1(1:neval) + 1i*varphi1(neval+1:neval2);
                % varphidot1_fast = [imag(HB_beta1)*real(psi1) + real(HB_beta1)*imag(psi1); -real(HB_beta1)*real(psi1) + imag(HB_beta1)*imag(psi1)]; 
                % Phi1 = [varphi1;varphidot1_fast];

                psi1 = topsi(Phi1);
                Phi1 = toPhi(beta1.*HB,psi1);
        
                if beta1 ~= 0
                    t = [0 1];
                    [T2,Phi] = ode15s(@(t,Phi)dPhi(t,Phi,XB_beta1),t,Phi1);
                    Phi_t = [Phi_t; Phi];
                    Phi2 = Phi(end, :).';
                else
                    Phi2 = Phi1;
                end

                % varphi2 = Phi2(1:neval2);
                % psi2 = varphi2(1:neval) + 1i*varphi2(neval+1:neval2);
                % varphidot2_fast = [imag(HP_gamma2)*real(psi2) + real(HP_gamma2)*imag(psi2); -real(HP_gamma2)*real(psi2) + imag(HP_gamma2)*imag(psi2)]; 
                % Phi2 = [varphi2;varphidot2_fast];

                psi2 = topsi(Phi2);
                Phi2 = toPhi(gamma2.*HP,psi2);


                if gamma2 ~= 0
                    t = [0 1];
                    [T3,Phi] = ode15s(@(t,Phi)dPhi(t,Phi,XP_gamma2),t,Phi2); %\varphi
                    Phi_t = [Phi_t; Phi];
                    Phi3 = Phi(end, :).';
                else
                    Phi3 = Phi2;
                end
                
                % varphi3 = Phi3(1:neval2);
                % psi3 = varphi3(1:neval) + 1i*varphi3(neval+1:neval2);
                % varphidot3_fast = [imag(HB_beta2)*real(psi3) + real(HB_beta2)*imag(psi3); -real(HB_beta2)*real(psi3) + imag(HB_beta2)*imag(psi3)]; 
                % Phi3 = [varphi3;varphidot3_fast];

                psi3 = topsi(Phi3);
                Phi3 = toPhi(beta2.*HB,psi3);

                if beta2 ~= 0
                    t = [0 1];
                    [T4,Phi] = ode15s(@(t,Phi)dPhi(t,Phi,XB_beta2),t,Phi3);
                    Phi_t = [Phi_t; Phi];
                    Phi4 = Phi(end, :).';
                else
                    Phi4 = Phi3;
                end

                % varphi4 = Phi4(1:neval2);
                % psi4 = varphi4(1:neval) + 1i*varphi4(neval+1:neval2);

                psi4 = topsi(Phi4);
               
                obj = psi4'*HP*psi4;
                objs(beta2idx,gamma2idx,beta1idx,gamma1idx) = obj;
            end
        end
    end
end

%[varphi1; varphi2; varphi3; varphi4]

psi4'*HP*psi4

% psi_t = []
obj_t = [];
for i = 1:size(Phi_t,1)
    Phi_i = Phi_t(i,1:neval2).';
    varphi_i = Phi_i(1:neval2);
    psi_i = varphi_i(1:neval) + 1i*varphi_i(neval+1:neval2);
    obj_i = psi_i'*HP*psi_i;
    obj_t = [obj_t obj_i];
end


T1_seq = T1 + 0;
T2_seq = T2 + T1_seq(end);
T3_seq = T3 + T2_seq(end);
T4_seq = T4 + T3_seq(end);

%[T1; T2+T1(end); T3+T2(end); T4+T3(end)]
%size([T1_seq;T2_seq;T3_seq;T4_seq])


[X,Y,Z,J] = ndgrid(beta2s,gamma2s,beta1s,gamma1s);
fSAVE = [X(:) Y(:) Z(:) J(:) objs(:)];
size(fSAVE)

txt1 = sprintf('maxcut_2ndorder_p2_8_qubit_new.txt');
fid1 = fopen(txt1,'w');
fprintf(fid1,'%1f %1f %1f %1f %1f\n', fSAVE.');
fclose(fid1);

figure(5)
plot(T1_seq, Phi_t(1:length(T1_seq),1), 'LineWidth',2)
hold on
plot(T2_seq, Phi_t(length(T1_seq)+1:length([T1_seq;T2_seq]),1), 'LineWidth',2)
hold on
plot(T3_seq, Phi_t(length([T1_seq;T2_seq])+1:length([T1_seq;T2_seq;T3_seq]),1), 'LineWidth',2)
hold on
plot(T4_seq, Phi_t(length([T1_seq;T2_seq;T3_seq])+1:length([T1_seq;T2_seq;T3_seq;T4_seq]),1), 'LineWidth',2)
% hold on
%plot([T1_seq;T2_seq;T3_seq;T4_seq], Phi_t(:,1), 'LineWidth',2);

xlabel('$t$','Interpreter','latex','FontSize',25)
ylabel('$C_i$','Interpreter','latex','FontSize',25)
set(gca,'FontSize',20)

ax = gca;
ax.FontSize = 27; %% tick font size
ax.FontName ='Times New Roman';
legend('Re$(C_0)$','Re$(C_1)$','Im$(C_0)$','Im$(C_1)$','location', 'best','Interpreter','latex','FontSize',10)
title('Probability Amplitudes','Interpreter','latex')
% print -dpdf probability_amplitudes



figure(6)
plot(T1_seq, Phi_t(1:length(T1_seq),1000), 'LineWidth',2)
hold on
plot(T2_seq, Phi_t(length(T1_seq)+1:length([T1_seq;T2_seq]),1000), 'LineWidth',2)
hold on
plot(T3_seq, Phi_t(length([T1_seq;T2_seq])+1:length([T1_seq;T2_seq;T3_seq]),1000), 'LineWidth',2)
hold on
plot(T4_seq, Phi_t(length([T1_seq;T2_seq;T3_seq])+1:length([T1_seq;T2_seq;T3_seq;T4_seq]),1000), 'LineWidth',2)
% hold on
% plot([T1_seq;T2_seq;T3_seq;T4_seq], Phi_t(:, 1000), 'LineWidth',2);

xlabel('$t$','Interpreter','latex','FontSize',25)
ylabel('$C_i$','Interpreter','latex','FontSize',25)
set(gca,'FontSize',20)

ax = gca;
ax.FontSize = 27; %% tick font size
ax.FontName ='Times New Roman';
legend('Re$(C_0)$','Re$(C_1)$','Im$(C_0)$','Im$(C_1)$','location', 'best','Interpreter','latex','FontSize',10)
title('Probability Amplitudes','Interpreter','latex')


figure(7)
plot(T1_seq, obj_t(1,1:length(T1_seq)), 'LineWidth',2)
hold on
plot(T2_seq, obj_t(1,length(T1_seq)+1:length([T1_seq;T2_seq])), 'LineWidth',2)
hold on
plot(T3_seq, obj_t(1,length([T1_seq;T2_seq])+1:length([T1_seq;T2_seq;T3_seq])), 'LineWidth',2)
hold on
plot(T4_seq, obj_t(1,length([T1_seq;T2_seq;T3_seq])+1:length([T1_seq;T2_seq;T3_seq;T4_seq])), 'LineWidth',2)
% hold on
% plot([T1_seq;T2_seq;T3_seq;T4_seq], obj_t, 'LineWidth',2);

xlabel('$t$','Interpreter','latex','FontSize',25)
ylabel('$C_i$','Interpreter','latex','FontSize',25)
set(gca,'FontSize',20)

ax = gca;
ax.FontSize = 27; %% tick font size
ax.FontName ='Times New Roman';
legend('Re$(C_0)$','Re$(C_1)$','Im$(C_0)$','Im$(C_1)$','location', 'best','Interpreter','latex','FontSize',10)
title('Probability Amplitudes','Interpreter','latex')



    
% function rhs = dPhi(t,Phi,hbar,X)
%     rhs = -X*Phi;

function con = check_L0(M)
    con_1 = full(all(diag(M)>=0));
    con_2 = full(all(M(~eye(size(M))<=0)));
    con = con_1*con_2;

function H = mat_toreal(H)
    H = [real(H) -imag(H);imag(H) real(H)];

function v = vec_toreal(v)
    v = [real(v);imag(v)];