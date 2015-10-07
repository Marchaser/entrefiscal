% collapse to state only
bp = 4*EpsilonPts*ZetaPts*APts-1;
yidx = [1:3*EpsilonPts*ZetaPts*APts bp+1:bp+8];
OmegaIdx = [3*EpsilonPts*ZetaPts*APts+1:4*EpsilonPts*ZetaPts*APts-1 bp+9:bp+12];

D = C(yidx);
F = G1(yidx,OmegaIdx);

C1 = C(OmegaIdx);
A1 = G1(OmegaIdx,OmegaIdx);
A2 = G1(OmegaIdx,yidx);

C2 = C1+A2*D;
A = A1+A2*F;

% eigen vector of A C2
SsNew = (eye(size(A,1))-A)\C2;
save('sys_helper','sys_helper');
SsOmega = Ss(OmegaIdx);
