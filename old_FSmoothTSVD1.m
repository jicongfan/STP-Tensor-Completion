function [X A B]=old_FSmoothTSVD1(Y,Omega,R,beta,lambda,eta,r,opt)

maxIter = opt.maxIter;

max_mu = opt.max_mu;

rho = opt.rho;

mu = opt.mu;

tol = opt.tol;

[n1,n2,n3] = size(Y);

Mask = zeros(n1,n2,n3);

Mask(Omega) = 1;

X = zeros(n1,n2,n3);

X(Omega) = Y(Omega);

X(~logical(Mask)) = mean(Y(logical(Mask)));

 

%initlize A B M(i) N(i) A_i B_i LAi LBi E

% XF = fft(X,[],3);

% AF = zeros(n1,R,n3);

% BF = zeros(n2,R,n3);

% for  i =1:n3

%     [UF SF VF]=svd(XF(:,:,i),'econ');

%     AF(:,:,i) = UF(:,1:R)*diag(sqrt(diag(SF(1:R,1:R))));

%     BF(:,:,i) = VF(:,1:R)*diag(sqrt(diag(SF(1:R,1:R))));

% end

 

%A = ifft(AF,[],3);

A =randn(n1,R,n3);

A1 = A;

A2 = A;

A3 = A;

 

AA = A;

 

%B = ifft(BF,[],3);

B =randn(n2,R,n3);

B1 = B;

B2 = B;

B3 = B;

 

BB = B;

 

E =zeros(n1,n2,n3);

 

LA1 = DiffM(size(A,1));

LA2 = DiffM(size(A,2));

LA3 = DiffM(size(A,3));

 

LB1 = DiffM(size(B,1));

LB2 = DiffM(size(B,2));

LB3 = DiffM(size(B,3));

 

MA1 = zeros(size(A,1),size(A,2)*size(A,3));

MA2 = zeros(size(A,2),size(A,1)*size(A,3));

MA3 = zeros(size(A,3),size(A,1)*size(A,2));

 

MB1 = zeros(size(B,1),size(B,2)*size(B,3));

MB2 = zeros(size(B,2),size(B,1)*size(B,3));

MB3 = zeros(size(B,3),size(B,1)*size(B,2));

 

% multipiler T PAi PBi QAi QBi

 

T = zeros(size(Y));

 

PA1 = zeros(size(A1,1),size(A1,2)*size(A1,3));

PA2 = zeros(size(A2,2),size(A2,1)*size(A2,3));

PA3 = zeros(size(A3,3),size(A3,1)*size(A3,2));

 

PB1 = zeros(size(B1,1),size(B1,2)*size(B1,3));

PB2 = zeros(size(B2,2),size(B2,1)*size(B2,3));

PB3 = zeros(size(B3,3),size(B3,1)*size(B3,2));

 

QA1 = zeros(size(A1));

QA2 = zeros(size(A2));

QA3 = zeros(size(A3));

 

 

QB1 = zeros(size(B1));

QB2 = zeros(size(B2));

QB3 = zeros(size(B3));

 

WA = zeros(size(AA));

WB = zeros(size(BB));

 

iter = 0;

 

while iter < maxIter

    

    iter = iter + 1;

    %update MAi MBi

    MA1 = prox_l1((mu*(LA1*double(tenmat(A1,1)))+PA1)/(mu+beta(1)),lambda(1)/(mu+beta(1)));

    MA2 = prox_l1((mu*(LA2*double(tenmat(A2,2)))+PA2)/(mu+beta(2)),lambda(2)/(mu+beta(2)));

    MA3 = prox_l1((mu*(LA3*double(tenmat(A3,3)))+PA3)/(mu+beta(3)),lambda(3)/(mu+beta(3)));

 

    MB1 = prox_l1((mu*(LB1*double(tenmat(B1,1)))+PB1)/(mu+beta(1)),lambda(1)/(mu+beta(1)));

    MB2 = prox_l1((mu*(LB2*double(tenmat(B2,2)))+PB2)/(mu+beta(2)),lambda(2)/(mu+beta(2)));

    MB3 = prox_l1((mu*(LB3*double(tenmat(B3,3)))+PB3)/(mu+beta(3)),lambda(3)/(mu+beta(3)));

 

      

%       MA1 = ClosedWL1(((LA1*double(tenmat(A1,1)))+PA1/mu),lambda(1)/mu,eps);

%       MA2 = ClosedWL1(((LA2*double(tenmat(A2,2)))+PA2/mu),lambda(2)/mu,eps);

%       MA3 = ClosedWL1(((LA3*double(tenmat(A3,3)))+PA3/mu),lambda(3)/mu,eps);

% 

%       MB1 = ClosedWL1(((LB1*double(tenmat(B1,1)))+PB1/mu),lambda(1)/mu,eps);

%       MB2 = ClosedWL1(((LB2*double(tenmat(B2,2)))+PB2/mu),lambda(2)/mu,eps);

%       MB3 = ClosedWL1(((LB3*double(tenmat(B3,3)))+PB3/mu),lambda(3)/mu,eps);

 

    

    % update Ai Bi

    temp_A1 =(inv(LA1'*LA1+eye(size(A1,1)))*((LA1'*(MA1-PA1/mu))+double(tenmat(A-QA1/mu,1))));

     for k = 1:size(A1,1)

         A1(k,:,:) = reshape(temp_A1(k,:),size(A1,2),size(A1,3));

     end

     

    temp_A2 =(inv(LA2'*LA2+eye(size(A2,2)))*((LA2'*(MA2-PA2/mu))+double(tenmat(A-QA2/mu,2))));

    for k = 1:size(A2,2)

         A2(:,k,:) = reshape(temp_A2(k,:),size(A2,1),size(A2,3));

    end

     

    temp_A3 =(inv(LA3'*LA3+eye(size(A3,3)))*((LA3'*(MA3-PA3/mu))+double(tenmat(A-QA3/mu,3))));

    for k = 1:size(A3,3)

         A3(:,:,k) = reshape(temp_A3(k,:),size(A3,1),size(A3,2));

    end

    

    temp_B1 =(inv(LB1'*LB1+eye(size(B1,1)))*((LB1'*(MB1-PB1/mu))+double(tenmat(B-QB1/mu,1))));

    for k = 1:size(B1,1)

         B1(k,:,:) = reshape(temp_B1(k,:),size(B1,2),size(B1,3));

     end

     

    temp_B2 =(inv(LB2'*LB2+eye(size(B2,2)))*((LB2'*(MB2-PB2/mu))+double(tenmat(B-QB2/mu,2))));

    for k = 1:size(B2,2)

         B2(:,k,:) = reshape(temp_B2(k,:),size(B2,1),size(B2,3));

    end

    

    temp_B3 =(inv(LB3'*LB3+eye(size(B3,3)))*((LB3'*(MB3-PB3/mu))+double(tenmat(B-QB3/mu,3))));

    for k = 1:size(B3,3)

         B3(:,:,k) = reshape(temp_B3(k,:),size(B3,1),size(B3,2));

    end

    

    % update AA  BB

    AA = reshape(prox_nuclear(double(tenmat(A-WA/mu,1)),eta/mu),n1,R,n3);

    BB = reshape(prox_nuclear(double(tenmat(B-WB/mu,1)),eta/mu),n2,R,n3);

 

    

    % update A  B

    A = tprod((AA+A1+A2+A3+(QA1+QA2+QA3+WA)/mu+tprod(Y.*Mask-E+T/mu,B)),tinv(tprod(tran(B),B)+4*teye(R,n3)));

    B = tprod(BB+B1+B2+B3+(QB1+QB2+QB3+WB)/mu+tprod(tran(Y.*Mask-E+T/mu),A),tinv(tprod(tran(A),A)+4*teye(R,n3)));

    

    % update E;

    E = Y.*Mask - tprod(A,tran(B))+T/mu;

    E(Omega) = 0 ;

    

    % update Multiplier

     leq = Y.*Mask - tprod(A,tran(B)) -E;

     

     leqA1 = (LA1*double(tenmat(A1,1)))-MA1;

     leqA2 = (LA2*double(tenmat(A2,2)))-MA2;

     leqA3 = (LA3*double(tenmat(A3,3)))-MA3;

     

     leqB1 = (LB1*double(tenmat(B1,1)))-MB1;

     leqB2 = (LB2*double(tenmat(B2,2)))-MB2;

     leqB3 = (LB3*double(tenmat(B3,3)))-MB3;

     

     geqA1 = A1 - A;

     geqA2 = A2 - A;

     geqA3 = A3 - A;

     

     geqB1 = B1 - B;

     geqB2 = B2 - B;

     geqB3 = B3 - B;

     

     geqAA = AA - A;

     geqBB = BB - B;

     

     

 

     err1 = max(abs(leq(:)));

     

     errDA = max([max(abs(leqA1(:))),max(abs(leqA2(:))),max(abs(leqA3(:)))]);

     errDB = max([max(abs(leqB1(:))),max(abs(leqB2(:))),max(abs(leqB3(:)))]);

 

     

     errA = max([max(abs(geqA1(:))),max(abs(geqA2(:))),max(abs(geqA3(:)))]);

     errB = max([max(abs(geqB1(:))),max(abs(geqB2(:))),max(abs(geqB3(:)))]);

 

     errAA = max(geqAA(:));

     errBB = max(geqBB(:));

     

     err=[err1,errDA,errDB,errA,errB,errAA,errBB];

    

    fprintf('    Iter%d\n   ',iter);

    fprintf('    Y-A*B-E %7.10f\n   ',err1);

    fprintf('    TVAi %7.10f\n   ',errDA);

    fprintf('    TVBi %7.10f\n   ',errDB);

    fprintf('    Ai %7.10f\n   ',errA);

    fprintf('    Bi %7.10f\n   ',errB);

    fprintf('    AA %7.10f\n   ',errAA);

    fprintf('    BB %7.10f\n   ',errBB);

    

    if max(err) < tol

        break;

    else

     T = T + mu*leq;

     PA1 = PA1 + mu*leqA1;

     PA2 = PA2 + mu*leqA2;

     PA3 = PA3 + mu*leqA3;

     

     PB1 = PB1 + mu*leqB1;

     PB2 = PB2 + mu*leqB2;

     PB3 = PB3 + mu*leqB3;

     

     QA1 = QA1 + mu*geqA1;

     QA2 = QA2 + mu*geqA2;

     QA3 = QA3 + mu*geqA3;

 

     QB1 = QB1 + mu*geqB1;

     QB2 = QB2 + mu*geqB2;

     QB3 = QB3 + mu*geqB3;

     

     WA = WA + mu*geqAA;

     WB = WB + mu*geqBB;

     

     

     mu = min(max_mu,mu*rho); 

    end

    

end

X = tprod(A,tran(B));

X(Omega) = Y(Omega);

 