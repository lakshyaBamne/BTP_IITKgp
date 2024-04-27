% LLF = load('Slowly_Moving_Shock_LLF_200.txt','-ascii');
CIDRA = load('Slowly_Moving_Contant_CIDRA_200_C2.txt','-ascii');
CIDRA_Old = load('Slowly_Moving_Contant_CIDRA_200_OLD.txt','-ascii');
% x = load('Entropy_CIDRA_400.txt','-ascii');
figure,
plot(CIDRA(:,1),CIDRA(:,2),'k-o')
hold on
plot(CIDRA_Old(:,1), CIDRA_Old(:,2),'r-s')
% xlim([-5,5]);
legend('CIDRA','OLD')