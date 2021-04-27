function intra = intra_classe (data,vecteur_centre,valeur_centre,nb_centre)
j=1; 
     for j=1:nb_centre
               cluster = find (vecteur_centre == j);
               for h=1:length(cluster)
                   for g=1:2
                       V(h,g) = data(cluster(h),g);
                   end
               end 
      d = DistMatrix (V,valeur_centre);    
      dict_min =  min(d,[],2);
      intra(j) = mean (dict_min);
      j = j+1;     
     end
end