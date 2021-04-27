clc ; 
clear ;

load yeastdata.mat

% pour réduire la taille de l'ensemble de données à 
% un sous-ensemble contenant uniquement les gènes les plus significatifs.

emptySpots = strcmp('EMPTY',genes);
yeastvalues(emptySpots,:) = [];
genes(emptySpots) = [];
 

% La fonction "isnanest" est utilisée pour identifier les gènes qui ont des données manquantes et
% des commandes d'indexation sont utilisées pour supprimer les gènes avec des données manquantes.

nanIndices = any(isnan(yeastvalues),2);
yeastvalues(nanIndices,:) = [];
genes(nanIndices) = [];


% " genevarfilter" est une fonction pour filtrer les gènes qui ont une petite variance dans le temps

mask = genevarfilter(yeastvalues);
yeastvalues = yeastvalues(mask,:);
genes = genes(mask);
numel(genes)

% clustering 
nb_centre = 2 ;

[vecteur_centre,valeur_centre] = kmeans(yeastvalues,nb_centre);
inter = inter_classe (valeur_centre);
eva1 = evalclusters(yeastvalues,vecteur_centre, 'DaviesBouldin');
eva2 = evalclusters(yeastvalues,vecteur_centre, 'DaviesBouldin');
% eva1 = evalclusters(yeastvalues,vecteur_centre, 'CalinskiHarabasz');
% eva2 = evalclusters(yeastvalues,vecteur_centre, 'CalinskiHarabasz');
% eva1 = evalclusters(yeastvalues,vecteur_centre, 'silhouette');
% eva2 = evalclusters(yeastvalues,vecteur_centre, 'silhouette');
DB = eva1.CriterionValues ;

% hybridation k-means et GLS

nbr_iteration = 0 ;
y = 1 ;
p = zeros(50) ;
q = 0 ;

while nbr_iteration < 10
     while 1
       affichage(nb_centre-1,1) = DB ;
       affichage(nb_centre-1,2) = nb_centre ;
       
            if eva2.CriterionValues > DB 
%               if eva2.CriterionValues < DB 
                for m=1:length(inter)
                    for s=1:length(inter)
                        penalisation(m,s) =  inter(m,s) *p(m,s);
                    end
                end 
                  
                     DB = DB + (y*sum(sum(penalisation)));  
%                        DB = DB - (y*sum(sum(penalisation)));  
 
          end
          
        nb_centre = nb_centre + 1 ;
        [vecteur_centre,valeur_centre] = kmeans(yeastvalues,nb_centre);
        inter2 = inter_classe (valeur_centre');
          eva2 = evalclusters(yeastvalues,vecteur_centre,'DaviesBouldin'); 
%           eva2 = evalclusters(yeastvalues,vecteur_centre,'CalinskiHarabasz');
%           eva2 = evalclusters(yeastvalues,vecteur_centre,'silhouette');
                           

                                   if eva2.CriterionValues < DB
                                        DB = eva2.CriterionValues ;
                                        inter = inter2 ;
                                   else
                                       q = q +1; 
                                       minima_locaux(nb_centre-1) = DB ;
                                       break ;
                                   end
     end
       
          for i=1:length(inter)
                for j=1:length(inter)
                    util(i,j) = inter(i,j)/ (1+ p(i,j));
                end 
          end 

[k,Indce] = max(inter(:));
[Indce_ligne, Indce_colone] = ind2sub(size(inter),Indce);
p(Indce_ligne,Indce_colone) = p(Indce_ligne,Indce_colone)+1 ;
nbr_iteration = nbr_iteration + 1;
 

end

 minn = 100; 
 for i=1:length(minima_locaux)
     if minn > (minima_locaux(i))&& (minima_locaux(i)>0)
         minn = minima_locaux(i);
         opt = i ;
     end
 end
 minn
 opt

% L'affichage 

figure
plot( affichage(:,2), affichage(:,1),'r-o');
title('Resultat de DaviesBouldin avec GLS & KHM');
xlabel('nombres  des  centres');
ylabel('Davies Bouldin');


