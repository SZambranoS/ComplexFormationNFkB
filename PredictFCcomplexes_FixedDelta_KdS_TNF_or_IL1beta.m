%This script allows to predict the number of activating complex based on
%transcriptomic data. For further details and citations, refer to
%Kizilirmak et al. https://www.biorxiv.org/content/10.1101/2021.12.07.471485v2


clear

hold off

close all


%RNA-seq data provided in the biorxiv manuscript 

[A,B,C]=xlsread('log.norm.rpkm.xls');


indexGTNF=[22:2:30];
indexKTNF=[12:2:20];
indexFTNF=[32:2:40];



indexGUT=[21:2:29];
indexRUT=[11:2:19];
indexBUT=[31:2:39];


%To estimate the number of complexes upon TNF, uncomment the following code
CuratedListActivatorsTNF; %List of proteins in the activating complex
indexG=indexGUT;
indexR=indexRUT;
indexB=indexBUT;
titlestr=('TNF-\alpha activators')
N=length(GeneID)
%This is the difference in order of magnitude of protein copies between TNF receptor and the
%remaining proteins in the complex. 
Delta=3
OoMin=[2,(2+Delta)*ones(1,N-1)]; %Minimum order of magnitude of proteins abundances in the list, from element 1 to N
OoMax=[3,(3+Delta)*ones(1,N-1)]; %Maximum order of magnitude of protein abunancess in the list, from element 1 to N



%%To estimate the number of complexes upon IL1-beta, uncomment the
%%following code
% CuratedListActivatorsIL1beta; %List of proteins in the activating complex
% indexG=indexGUT;
% indexK=indexRUT;
% indexF=indexBUT;% indexG=indexGUT;
% titlestr=('IL-1\beta activators')
% N=length(GeneID)
% Delta=2
% OoMin=[3, 2, (3+Delta)*ones(1,N-2)]; %Minimum order of magnitude of proteins abundances in the list, from element 1 to N
% OoMax=[4, 3, (4+Delta)*ones(1,N-2)];  %Maximum order of magnitude of protein abunancess in the list, from element 1 to N
% 


%Here we define how much we allow to vary the dissociation constant K_D of
%each biochemical reaction. 
%Maximum orders of magnitudes of the KD considered

minOoMKD=0 %minimum order of magnitude considered for K_D.

maxOoMKD=3 %Maximum order of magnitude considered for K_D.

Nsims=1000; %Number of simulations performed



%We find the transcriptional differences.

for n=1:length(GeneID)
    
    indexuseful=find(strcmp(B(:,1),GeneID{n})); %FIND THE GENE
    indexgene=indexuseful-1;
    GeneValues=2.^A(indexgene,:);
    
    Gvalues=GeneValues(indexGUT);
    Rvalues=GeneValues(indexRUT);
    Bvalues=GeneValues(indexBUT);
    
    
    x=[Bvalues',Rvalues',Gvalues'];
    rnavalues=[Bvalues,Rvalues,Gvalues]
  
    
    Tablenorm(n,:)=rnavalues/mean(rnavalues);
    
    
    
end;


AllratiosBtoRsim=[];
AllratiosGtoRsim=[];


for nsims=1:Nsims
    
    
    
    Tablesimulatedprot=[];
    
    for nprot=1:length(OoMin)
        %This is how much we allow to vary the number of protein copies of
        %each element of the list, between the orders of magnitudes
        %selected
        vectorproteinfactor(nprot)=10^(OoMin(nprot)+(OoMax(nprot)-OoMin(nprot))*rand);
        Tablesimulatedprot(nprot,:)=Tablenorm(nprot,:)*vectorproteinfactor(nprot);
        
    end;
    
    
    for nprot=1:length(OoMin)-1
        %For each biochemical reaction we select a K_D in the range
        %prescribed. 
        
        listkD(nprot)=10^(minOoMKD+(maxOoMKD-minOoMKD)*rand);
        
    end;
    
    ResultsKd=[];
    
    [nproteins,nreplicates]=size(Tablesimulatedprot);
    
    %This is our calculation of the equilibrium concentration.
    for nr=1:nreplicates
        listproteins=Tablesimulatedprot(:,nr);
        ResultsKd(nr)=sequentialassociationreaction(listproteins,listkD);
        
    end;
    
    
    BUT=ResultsKd(1:5); %Equilibrium derived from the replicates of B
    RUT=ResultsKd(6:10); %Equilibrium derived from the replicates of R
    GUT=ResultsKd(11:15); %Equilibrium derived from the replicates of G
    
    m=1;
    
    ratioBtoRsim=[];
    ratioGtoRsim=[];
    
    
    %We calculate all the possible ratios B vs R and G vs R
    
    for nrep=1:length(BUT)
        
        
        for k=1:length(GUT)
            
            ratioBtoRsim(m)=BUT(nrep)/RUT(k);
            ratioGtoRsim(m)=GUT(nrep)/RUT(k);
            
            m=m+1;
            
            
        end
        
        
    end;
    
    
    AllratiosBtoRsim=[AllratiosBtoRsim,ratioBtoRsim];
    AllratiosGtoRsim=[AllratiosGtoRsim,ratioGtoRsim];
    
    
    
    
end;


histogram(AllratiosBtoRsim,'facecolor','b')
hold on
histogram(AllratiosGtoRsim,'facecolor','g')
