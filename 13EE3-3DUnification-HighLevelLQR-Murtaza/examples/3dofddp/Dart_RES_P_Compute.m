clear;
close all;
%%
F1 = [zeros(3,3) eye(3);zeros(3,3) zeros(3,3)];
F2 = [zeros(3,3) eye(3);zeros(3,3) zeros(3,3)];
F3 = [zeros(3,3) eye(3);zeros(3,3) zeros(3,3)]; 
F4 = [zeros(3,3) eye(3);zeros(3,3) zeros(3,3)];
F5 = [zeros(3,3) eye(3);zeros(3,3) zeros(3,3)];
% F5_Full = [0 1;0 0];
% F6 = [zeros(2,2) eye(2);zeros(2,2) zeros(2,2)];
% F6 = zeros(1);
% F7 = zeros(1);
% F8 = zeros(1);

%%
% F_Speed_Full = blkdiag(F1,F2,F3,F4,F5,F6,F7,F8);
F_Speed_Full = blkdiag(F1,F2,F3,F4,F5);
%%
G1 = [zeros(3,3);eye(3,3)];
G2 = [zeros(3,3);eye(3,3)];
G3 = [zeros(3,3);eye(3,3)];
G4 = [zeros(3,3);eye(3,3)];
G5 = [zeros(3,3);eye(3,3)];
G6 = eye(1);
G7 = eye(1);
G8 = eye(1);

%%
% G_Speed_Full = blkdiag(G1,G2,G3,G4,G5,G6,G7,G8);
G_Speed_Full = blkdiag(G1,G2,G3,G4,G5);
 %%
 KpEE = 750;
 KvEE = 250;
 
 KpOr = 150.0;
 KvOr = 50;
 
 KpCOM = 5;
 KvCOM = 2;
 KvSpeedReg = 2;
 
 KpEEL = eye(3)*KpEE;
 KvEEL = eye(3)*KvEE;
 
 KpEER = eye(3)*KpEE;
 KvEER = eye(3)*KvEE;
 
 
 KpOrL = eye(3)*KpOr;
 KvOrL = eye(3)*KvOr;
 
 KpOrR = eye(3)*KpOr;
 KvOrR = eye(3)*KvOr;
 
 KpCOM1 = eye(3)*KpCOM;
 KvCOM1 = eye(3)*KvCOM;
 
%  [-KpTh -KvTh],[-KpEE -KvEE],[-KpOr -KvOr],[-KvReg*eye(2)]
 K =  blkdiag([-KpEEL -KvEEL], [-KpEER -KpEER],[-KpOrL -KvOrL],[-KpOrR -KvOrR],[-KpCOM -KvCOM],[-KvSpeedReg*eye(3)]);  
% K =  blkdiag([-KpEEL -KvEEL], [-KpEER -KpEER],[-KpOrL -KvOrL],[-KpOrR -KvOrR],[-KpTh -KvTh],[-KvSpeedReg*eye(5)]);  
 
%  Fcl = F_Speed_Full+G_Speed_Full*K;
 
 %%
 
% QQ_Pose = eye(size(F_Speed_Full));
% QQ_Speed = eye(35);
QQ_Full = eye(size(F_Speed_Full));

% QQ = diag([10,10,1,1,1,1,1,1,0.1000,0.1000,0.1000,0.1000,0.1000,0.1000]);
% P1 = lyap(Fcl,QQ);

% eta = [th-thd;dth-dthd;rEE-rEEd;drEE-drEEd;Or-Ord;dOr-dOrd;dq(4:5)];
% eta = [th-thd;dth-dthd;rEE-rEEd;drEE-drEEd;Or-Ord;dOr-dOrd;dq];

[P_Speed_Full,K,~,INFO] = icare(F_Speed_Full,G_Speed_Full,QQ_Full);

% P_Speed_Full = lyap(Fcl,QQ_Full);
P1_Full = round(P_Speed_Full,5);

writematrix(P1_Full,'P_space_SpeedReg_Full.txt','Delimiter',' ')
writematrix(F_Speed_Full,'F_space_SpeedReg_Full.txt','Delimiter',' ')
writematrix(G_Speed_Full,'G_space_SpeedReg_Full.txt','Delimiter',' ')


% LfV_x = eta'*(F'*P+P*F)*eta;
% LgV_x = 2*eta'*P*G;
% V_x = eta'*P*eta;

% lambda_minQ = min(eig(QQ_Speed));
% lambda_maxP = max(eig(P1_Speed));

% lambda_minQ = min(eig(QQ_Pose));
% lambda_maxP = max(eig(P1_Pose));