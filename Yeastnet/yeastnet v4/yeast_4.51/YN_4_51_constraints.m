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
%NOTE: SINCE TCA OK WITH iND CONSTRAINTS, these may be uneeded?
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

%add H+ and H2O to biomass
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

%YAY! I now have YN aerobically producing biomass, using glycolysis and TCA
%cycles! THIS IS BETTER THAN iND DOES!... now, back to anaerobic headaches.
