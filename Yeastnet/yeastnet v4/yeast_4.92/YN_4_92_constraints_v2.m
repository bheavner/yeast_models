%BH 6/1/11

%script to constrain YN 4.92 to apply constraints from iMM904, then to
%check glycolysis, TCA, and ethanol production and biomass (anaerobic)

%known YN issues:
%1) the FBA solution does not have fluxes through the CoA biosynthetic 
% pathway (note: there is speculation about other CoA
% biosynthetic pathways in Yeast. doi:10.1016/j.plipres.2005.04.001 argues
% against such pathways. Specifically, CoA should not be synthesized in the
% mitochondrion.
%
%2) also related to CoA metabolism, YN is not a pantothenate auxotroph. 
% Per doi:10.1016/j.plipres.2005.04.001 (A VERY GOOD REFERENCE!):
% "Yeast is a pantothenate auxotroph, although the requirement for 
% pantothenate can be overcome by the addition of b-alanine to the growth 
% medium. The product of the FEN2 gene of Saccharomyces cerevisiae is a 
% hydrogen ion-pantothenate symporter belonging to the major facilitator 
% superfamily localized on the plasma membrane"
%
%3) anaerobic metabolism isn't clearly blocked b/c fatty acid desaturase 
% activty, as it should be. This is due to my bypassing "acyl-CoA", and
% unresolved issues with sterol biosynthesis and substitution.
%
%4) I'm suspicious of the phospholipid metabolism in the mitochondria
%
%5) phosphatidate is required to compensate for problems with the anaerobic
% biosynthesis of NADH and acyl-CoA.
%
%6) anaerobic FBA doesn't have all fluxes through glycolysis (but minimal
%1-norm FBA does). I take this to suggest that there remain yet more
%constraints that should be applied so that FBA better reflects internal
%fluxes - but that may also be due to structural issues with lipid, sterol, 
%and cofactor reconstruction.

%initcobratoolbox;
%YN=readcbmodel;
model=YN;
sln=optimizecbmodel(model) %0.0927
sln.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis not ok
sln.x([967 302 304 282 663 661 838 837 1028 1026 454 717]) %TCA ok
sln_min=optimizecbmodel(model,[],'one') %.0927
sln_min.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln_min.x([967 302 304 282 663 661 838 837 1028 1026 454 717])%TCA ok

rxns=[11 75 95 128 130 131 132 133 135 136 137 458 588 793 794 795 797 ...
    798 799 809 810 811 812 813 872 1135 1137 1158 1215 1280 1284];
%rxns is constraints from iMM904 with same genes as in YN - to to LB=0
%UB=1000. I took out rxn 1094 b/c it breaks biomass production

rxns2=[997 1000 1023 1024 1025]; %from iMM904, to LB=-1000 UB=0

rxns3=[408 736 920 1093]; %unconstrained in iMM904, but I think they should 
%be LB=0, UB=1000

model.lb(rxns)=0;
sln=optimizecbmodel(model) %0.0903
sln.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln.x([967 302 304 282 663 661 838 837 1028 1026 454 717]) %TCA ok
sln_min=optimizecbmodel(model,[],'one') %0
sln_min.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis not ok
sln_min.x([967 302 304 282 663 661 838 837 1028 1026 454 717])%TCA not ok

model.ub(rxns2)=0;
sln=optimizecbmodel(model) %0.0903
sln.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln.x([967 302 304 282 663 661 838 837 1028 1026 454 717]) %TCA ok
sln_min=optimizecbmodel(model,[],'one');
sln_min.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln_min.x([967 302 304 282 663 661 838 837 1028 1026 454 717])%TCA ok

model.lb(rxns3)=0; 
sln=optimizecbmodel(model) %0.0903
sln.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln.x([967 302 304 282 663 661 838 837 1028 1026 454 717]) %TCA ok
sln_min=optimizecbmodel(model,[],'one');
sln_min.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln_min.x([967 302 304 282 663 661 838 837 1028 1026 454 717])%TCA ok

%next, some rxn changes based on iMM904
%first, modify rxn 11 to be like iMM rxn 1270 (move proton)
model.S(799,11)=1; %was -1
sln=optimizecbmodel(model) %0.0903
sln.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln.x([967 302 304 282 663 661 838 837 1028 1026 454 717]) %TCA ok
sln_min=optimizecbmodel(model,[],'one');
sln_min.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln_min.x([967 302 304 282 663 661 838 837 1028 1026 454 717])%TCA ok

%next change rxns 75, 95, and 588 to use ATP instead of Pi (like iMM 1293
%and 1295)
model.S(1321,75)=0;model.S(438,75)=-1;model.S(394,75)=1; model.lb(75)=0;
model.S(1321,95)=0;model.S(438,95)=-1;model.S(394,95)=1; model.lb(95)=0;
model.S(1321,588)=0;model.S(438,588)=-1;model.S(394,588)=1; model.lb(588)=0;
sln=optimizecbmodel(model) %0.0903
sln.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln.x([967 302 304 282 663 661 838 837 1028 1026 454 717]) %TCA ok
sln_min=optimizecbmodel(model,[],'one');
sln_min.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln_min.x([967 302 304 282 663 661 838 837 1028 1026 454 717])%TCA ok

%next combine YN 1135 and 1137 like iMM 358
model.S([530 532],1135)=0;
model.S([530 796],1137)=1; model.S([532 791],1137)=-1;
model.lb(1137)=0;
sln=optimizecbmodel(model) %0.0903
sln.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln.x([967 302 304 282 663 661 838 837 1028 1026 454 717]) %TCA ok
sln_min=optimizecbmodel(model,[],'one');
sln_min.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln_min.x([967 302 304 282 663 661 838 837 1028 1026 454 717])%TCA ok

%next combine YN 1280 and 1284 like iMM 1557
model.S([1548 1550],1280)=0;
model.S([1548 796],1284)=1; model.S([1550 791],1284)=-1;
model.lb(1284)=0;
sln=optimizecbmodel(model) %0.0903
sln.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln.x([967 302 304 282 663 661 838 837 1028 1026 454 717]) %TCA ok
sln_min=optimizecbmodel(model,[],'one');
sln_min.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln_min.x([967 302 304 282 663 661 838 837 1028 1026 454 717])%TCA ok

%next YN 1158 should be like iMM 694
model.S(711,1158)=1;model.S(713,1158)=-1;
model.lb(1158)=0;
sln=optimizecbmodel(model) %0.0903
sln.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln.x([967 302 304 282 663 661 838 837 1028 1026 454 717]) %TCA ok
sln_min=optimizecbmodel(model,[],'one');
sln_min.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln_min.x([967 302 304 282 663 661 838 837 1028 1026 454 717])%TCA ok

%finally, a few more constraints from metacyc:
forward = [153 210 404 412 414 461 569 959 962];
reverse = [52 58 348 561 562 746 747 748 749];
model.lb(forward)=0;
model.ub(reverse)=0;
sln=optimizecbmodel(model) %0.0903
sln.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln.x([967 302 304 282 663 661 838 837 1028 1026 454 717]) %TCA ok
sln_min=optimizecbmodel(model,[],'one');
sln_min.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
sln_min.x([967 302 304 282 663 661 838 837 1028 1026 454 717])%TCA ok

%------------------ ANAEROBIC BELOW

%And the following is less supported by the literature, to make anaerobic work:

%To produce anaerobic biomass, set the following five to have UB=1000:
%* lipid pseudoreaction [no 14-demethyllanosterol, no 
%   ergosta-5,7,22,24(28)-tetraen-3beta-ol] 2125
%* ergosterol exchange 1768
%* lanosterol exchange 1923
%* zymosterol exchange 2122
%* phosphatidate exchange 2027

model.ub(2010)=0; %oxygen constraint (was 1000)

ana_sln=optimizecbmodel(model) %5.1e-25
ana_sln.x(1772) %-2 - it is making ethanol.
ana_sln.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis not ok
ana_sln.x([967 302 304 282 663 661 838 837 1028 1026 454 717]) %TCA not ok
ana_sln_min=optimizecbmodel(model,[],'one');
ana_sln_min.x(1772) %0 - it is not making ethanol.
ana_sln_min.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
ana_sln_min.x([967 302 304 282 663 661 838 837 1028 1026 454 717])%TCA ok

model.ub([2125 1768 1923 2122 2027])=1000;

ana_sln=optimizecbmodel(model) %0.0218
ana_sln.x(1772) %-1.7464 - it is making ethanol.
ana_sln.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
ana_sln.x([967 302 304 282 663 661 838 837 1028 1026 454 717]) %no TCA - ok
ana_sln_min=optimizecbmodel(model,[],'one');
ana_sln_min.x(1772) %-1.7464 - it is not making ethanol.
ana_sln_min.x([1173 536 469 892 452 488 898 899 368 968] ) %glycolysis ok
ana_sln_min.x([967 302 304 282 663 661 838 837 1028 1026 454 717])%no TCA - ok

%---------------

%so some weirdness with anaerobic glycolysis...hrmmm

