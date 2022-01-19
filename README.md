# Apprentissage-Grande-Dimension

On considère une matrice de données X et un vecteur d’observations Y à expliquer. Les observations portent sur p variables, mesurées sur n individus. Les données fil rouge de ce repo seront des données d’apprentissage supervisée en génomique. 

En présence de nombreuses variables explicatives, on suppose généralement que peu d’entre elles sont pertinentes pour modéliser/prédire Y. Il existe quatre familles de méthodes permettant de contourner le fléau de grande dimension : 
  - tests multiples : utilisés en pré-traitement pour filtrer les variables
  - réduction de dimension : utilisés en pré-traitement pour réduire la dimension de l’espace des variables
  - choix de modèles : pour choisir le meilleur sous-modèle
  - régressions sous contraintes (ou pénalisées) : pour contraindre le nombre de paramètres dans le modèle
