within Steps.TPComponents;

package Correlations
  // enumeration of avaliable Nusselt Number correlation
  type NuCorrType = enumeration(Gnielinski, Ngo, Liao, Xin);
  
  // friction factor correlation
  type FFCorrType = enumeration(Gnielinski, Ngo, Xin); 
       
end Correlations;
