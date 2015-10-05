[EpsilonSubMesh,ASubMesh] = ndgrid(EpsilonGrid,AGrid);
VManager2 = reshape(VManager, [EpsilonPts,ZetaPts,APts]);
VManager2 = squeeze(VManager2(:,5,:));
VWorker2 = reshape(VWorker, [EpsilonPts,ZetaPts,APts]);
VWorker2 = squeeze(VWorker2(:,5,:));
VMax2 = reshape(VMax, [EpsilonPts,ZetaPts,APts]);
VMax2 = squeeze(VMax2(:,5,:));

plot(AGrid,ADist);

[EpsilonSubMesh,ASubMesh] = ndgrid(EpsilonGrid,AGrid);
NWorker2 = reshape(NWorker,[EpsilonPts,ZetaPts,APts]);
NWorker2 = squeeze(NWorker2(:,5,:));
mesh(EpsilonSubMesh,ASubMesh,NWorker2);

[EpsilonSubMesh,ASubMesh] = ndgrid(EpsilonGrid,AGrid);
Occ2 = reshape(OccPolicy,[EpsilonPts,ZetaPts,APts]);
Occ2 = squeeze(Occ2(:,end,:));
mesh(EpsilonSubMesh,ASubMesh,Occ2);

[EpsilonSubMesh,ASubMesh] = ndgrid(EpsilonGrid,AGrid);
Dist2 = squeeze(Dist(:,7,:));
mesh(EpsilonSubMesh,ASubMesh,Dist2);

figure(1);
mesh(EpsilonSubMesh,ASubMesh,VManager2);

figure(2);
mesh(EpsilonSubMesh,ASubMesh,VWorker2);

figure(3);
mesh(EpsilonSubMesh,ASubMesh,VMax2);