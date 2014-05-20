%BH 5/16/11

%script to constrain YN 4.51 for glycolysis, TCA, pentose phosphate
%(aerobic) and ethanol production and biomass (anaerobic)

%known YN issues:
%1) the FBA solution does not have fluxes through the CoA biosynthetic 
% pathway (including, but maybe not limited to rxns 758 1942 1943 1944 
% 1949 1947 1945 1946) (note: there is speculation about other CoA
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
%2) anaerobic metabolism isn't clearly blocked b/c fatty acid desaturase 
% activty, as it should be. This is due to my bypassing "acyl-CoA", and
% unresolved issues with sterol biosynthesis and substitution.
%
%3) I'm suspicious of the phospholipid metabolism in the mitochondria
%
%4) phosphatidate is required to compensate for problems with the anaerobic
% biosynthesis of NADH and acyl-CoA.


%initcobratoolbox;
%YN=readcbmodel;
model=YN;
sln1=optimizecbmodel(model) %0.1676
sln1.x([1918 542 475 801 457 495 806 807 372 878])  %glycolysis not ok
sln1.x([877 306 308 285 591 589 749 748 931 461 929 463 459 649])%TCA not ok
%10 0s

rxns=[238 267 268 269 454 502 503 505 1399 1647 1648 1649 1650 1942 ...
1943 1944 1945 2119 2120 2121 2122]; %iND constrainted rxns with no gene

old_rxns=[1635 1047 121 1915 125 149 1902 512 935 249 253 254 255 1636 ...
1637 1638 1645 1646 1647 1648 1649 1650 277 267 268 269 279 346 1914 226 ... 
1354 411 110 656 2123 462 1918 511 512 513 515 534 572 378 606 604 1921 ...
2041 2042 2043 2044 2045 2046 2047 2048 346 1903 2062 2063 2064 2065 2066 ...
2067 2068 2069 2070 2071 2052 2053 2054 2055 2056 2057 2058 2059 2060 2061 ...
2039 1931 716 790 1720 1721 1722 1723 1724 1725 1634 1322 1732 988]; 
%iND constrained rxns with genes

model.lb(rxns)=0;
sln=optimizecbmodel(model) %0.1676 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolysis ok
sln.x([877 306 308 285 591 589 749 748 931 461 929 463 459 649])%TCA not ok
%7 0s

model.lb(old_rxns)=0;
sln=optimizecbmodel(model) %0.1676 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolysis ok
sln.x([877 306 308 285 591 589 749 748 931 461 929 463 459 649])%TCA not ok
%8 0s

model.lb([369 947 986 1715 1716 1911 1920 1923])=0; 
%my ATP constraints for TCA
%NOTE: SINCE TCA OK WITH iND CONSTRAINTS, these may be uneeded? - but evidence
%from metacyc supports
sln=optimizecbmodel(model) %0.1676 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolysis ok
sln.x([877 306 308 285 591 589 749 748 931 461 929 463 459 649])%TCA not ok
% 7 0s

model.S(518,2125)=0;model.S(517,2125)=-1; %fix peroxisome CoA transport
sln=optimizecbmodel(model) %0.1676 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolysis ok
sln.x([877 306 308 285 591 589 749 748 931 461 929 463 459 649])%TCA not ok
% 8 0s

model.lb([649 328 648 898 250 544 805 1921])=0; %NAD cytosol constraints
sln=optimizecbmodel(model) %0.0948 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolysis ok
sln.x([877 306 308 285 591 589 749 748 931 461 929 463 459 649])%TCA not ok
% 7 0s

model.ub([1912 173])=0; %NAD cytosol constraints
sln=optimizecbmodel(model) %0.0948 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolysis ok
sln.x([877 306 308 285 591 589 749 748 931 461 929 463 459 649])%TCA not ok
% 6 0s

model.lb([649])=0; %NAD mito constraints
sln=optimizecbmodel(model) %0.0948 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolysis ok
sln.x([877 306 308 285 591 589 749 748 931 461 929 463 459 649])%TCA not ok
% 6 0s

%add H+ and H2O to biomass def
model.S(1440,2126)=-59.276;
model.S(765,2126)=58.7162;
sln=optimizecbmodel(model) %0.0948 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolysis ok
sln.x([877 306 308 285 591 589 749 748 931 461 929 463 459 649])%TCA not ok
% 7 0s

model = addReaction(model,'ATPmaint','s_0446 + s_1434 -> s_0400 + s_1207 + s_0764');
model.lb(2128)=0;
sln=optimizecbmodel(model) %0.0948 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolysis not ok
sln.x([877 306 308 285 591 589 749 748 931 461 929 463 459 649])%TCA not ok
% 7 0s

%ATP maintanence does bad things?


model.S(770,591)=1; %rxn 591 should produce H+:
sln=optimizecbmodel(model) %0.0947 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolysis ok? 
%(rxn 475 has negative flux)
sln.x([877 306 308 285 591 589 749 748 931 461 929 463 459 649])%TCA not ok
% 6 0s

model.lb([308 285])=0; %was -1000, change per metacyc
sln=optimizecbmodel(model) %0.0947 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolysis ok? 
%(rxn 475 has negative flux)
sln.x([877 306 308 285 591 589 749 748 931 929 459 649])%TCA not ok
% 6 0s

model.lb([1904])=0; %was 1000 per sgd 5/12
sln=optimizecbmodel(model) %0.0927 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolosis ok!
sln.x([877 306 308 285 591 589 749 748 931 929 459 649])%TCA improved!
% 2 0s

model.S([499 838 1440],286)=0; %remove reaction 286, which is in the 
%mitochondria (as rxn 285)
sln=optimizecbmodel(model) %0.0927 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolosis ok!
sln.x([877 306 308 285 591 589 749 748 931 929 459 649])%TCA ok!

model.lb(273)=0; %phosphatidate transport - constrained per metacyc 5/12
sln=optimizecbmodel(model) %0.0927 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolosis ok!
sln.x([877 306 308 285 591 589 749 748 931 929 459 649])%TCA ok!

model.lb([407 558])=0; %direction on CoA from metacyc 5/13
model.ub([567 568])=0; %direction on CoA from metacyc 5/13
sln=optimizecbmodel(model) %0.0927 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolosis ok!
sln.x([877 306 308 285 591 589 749 748 931 929 459 649])%TCA ok!

%------------------

%And the following is less supported by the literature, to make anaerobic work:

model.ub([1791 1898 1844])=1000; %allow sterols(were 0)
sln=optimizecbmodel(model) %0.0939 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolosis ok!
sln.x([877 306 308 285 591 589 749 748 931 929 459 649])%TCA not ok!
%rxn 589 has 0 flux, but there's a flux through 591. and it can be 591 or 589 - no prob.
%It's NAD or NADP as cofactors. So this is more evidence that there's a NAD 
%issue (and it's related to sterols somehow).

model=addexchangerxn(model,{'s_1216'}); %add phosphatidate [cytoplasm] exchange
sln=optimizecbmodel(model) %0.0945 - not overconstrained!
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolosis ok!
sln.x([877 306 308 285 591 589 749 748 931 929 459 649])%TCA ok!

%remove '14-demethyllanosterol [cytoplasm]' and 
%'ergosta-5,7,22,24(28)-tetraen-3beta-ol [cytoplasm]' from the biomass 
%definition:
model.S([124 633],1951)=0; 
sln=optimizecbmodel(model) %0.0945
sln.x(1792) %0 - it is not making ethanol.
sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolysis ok
sln.x([877 306 308 285 591 589 749 748 931 929 459 649])%TCA ok!

model.ub(1861)=0; %oxygen constraint (was 1000)

ana_sln=optimizecbmodel(model) %0.0230
ana_sln.x(1792) %-1.7321 - it is making ethanol.
ana_sln.x([1918 542 475 801 457 495 806 807 372 878])  %glycolysis ok
ana_sln.x([877 306 308 285 591 589 749 748 931 929 459 649])%TCA ok - no flux!

%---------------

%This makes me pretty happy!

