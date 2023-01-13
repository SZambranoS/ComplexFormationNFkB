function Cs = sequentialassociationreaction(listproteins,listkd)

%This function gives the number of complexes formed by the successive
%association of the molecules in numbers of proteins as in listproteins and
%with equilibrium perscribed as the kds in listkd so listkd(j) is the kd of
%the reaction between listroteins(j) and listproteins(j+1)




for n=1:length(listproteins)-1

    kD=listkd(n);
    
    if n==1;
        C0=listproteins(n);
    end; 
       
    A0=C0; 
    B0=listproteins(n+1);
    
    [As,Bs,Cs]=AnalyticalSolutionGeneralComplex(kD,1,1,A0,B0); %To try different stoichiometries, use other values other than 1, 1
    
    C0=Cs;
    
end; 


