# TD1 : Rappels sur le modèle linéaire

# Exercice 1 : Petits exemples

# Exercice 1.1 : jeu de données iris

# 1.

data('iris')
dim(iris) # (150, 5)
head(iris)

# 2.

model <- lm(Petal.Length ~ ., data = iris) # Y = Petal.Length et X = variables restantes
model # estimateurs associés
summary(model) # sigificativité des estimateurs

# Toutes les variables sont retenues comme significatives par le modèle (avec une p-value < 0.05). Cependant,
# la variable principale semble être la longueur de la sépale Sepal.Length. On remarque que la variable
# Species apparaît deux fois. Ceci est dû au fait qu’il s’agit d’une variable discrète à trois états. Il faut donc
# la modéliser à l’aide de deux coefficients dans le vecteur \beta : l’un des états (ici setosa) est choisi comme base,
# et la variable Species est remplacée par deux indicatrices (versicolor et virginica). Les coefficients liés à
# ces variables indiquent alors la différence de valeur moyenne quand on passe de setosa à l’état correspondant.
# En général, on ne prend pas les variables qualitatives comme explicatives dans un modèle linéaire gaussien.

# 3.

iris$versicolor <- iris$Species == 'versicolor'
modellogit <- glm(versicolor ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
                  family = binomial, data = iris) # modèle de régression logistique
summary(modellogit)

# La variable "Sepal.Width" caractérise le mieux l'espcèce "versicolor".

# Exercice 1.2 : jeu de données airquality

# 1.

data("airquality")
dim(airquality) # (153, 6)
head(airquality)

# 2.

modelair <- lm(Ozone ~ ., data = airquality) # Y = Ozone et X = variables restantes
modelair
summary(modelair)

# Les coefficients les plus significatifs sont un coefficient positif associé 
# à la température et un coefficient négatif associé au vent, ce qui correspond 
# bien à nos à priori (plus il fait chaud, plus le taux d’ozone est élevé, et 
# inversement pour le force du vent).

# 3.

boxplot(Ozone ~ Month, data = airquality)

summer <- airquality[airquality$Month < 8, ]
modelair2 <- lm(Ozone ~ ., data = summer) # Y = Ozone et X = variables restantes
modelair2
summary(modelair2)

# Si on réduit les données aux mois de Mai, Juin et Juillet, les conclusions 
# restent identiques que sur l’ensemble des données.

res <- lm(Ozone ~ Month, data = summer)
summary(res)

# Si on explique le taux d’ozone par la variable mois, on obtient un coefficient 
# très significatif et positif. Cette fluctuation par la présence ou non 
# des variables vent et température laisse penser que la variable mois est 
# sans doute très corrélée à l’une ou l’autre de ces variables, ce qui
# rend l’interprétation des coefficients difficiles.

cor.test(summer$Wind, summer$Month)

cor.test(summer$Temp, summer$Month)

# Les variables semblent corrélées.

# Exercice 2 : Problème de colinéarité

# Exercice 2.1 : jeu de données Prostate

# 1.

library(lasso2)
data("Prostate")
dim(Prostate) # (97, 9)
head(Prostate)

modellcavol <- lm(lcavol ~ ., data = Prostate)
summary(modellcavol)

# Trois variables (age, lcp et lpsa) semblent expliquer la variable lcavol 
# (p-valeurs significatives).

# 2.

library(corrplot)

cov <- cor(Prostate)
corrplot(cov)

# La variable lpsa est très corrélée à la variable réponse lcavol, ce qui peut 
# expliquer pourquoi elle a une influence significative. Les autres variables, 
# pgg45, gleason, lcp et svi semblent être corrélées. Cela peut poser un problème 
# lorsque l’on interprète les variables significatives dans le modèle linéaire.

# 3.

library(car)

VIF <- vif(modellcavol)
VIF

# VIF < 10 => les variables ne sont pas colinéaires.

# Exercice 2.2 : passage à la grande dimension

# 1.

# (on se donne p-1 variables très corrélées à X^1 et p-1 très corrélées à X^{p+1})

library(MASS)

CreateData <- function(p, n, rho){# p : nb variables, n : taille de l'échantillon
  
  # créer la matrice de covariance
  sigma1 <- matrix(rho,p,p)
  diag(sigma1) <- 1
  sigma2 <- rbind(cbind(sigma1, matrix(0,p,p)), cbind(matrix(0,p,p), sigma1)) # concaténation 2p+1
  
  # créer la matrice X
  X <- mvrnorm(n = n, mu = rep(0, nrow(sigma2)), sigma2)
  
  # Y = X^1 + X^{p+1} + epsilon
  data <- cbind(X, X[,1] + X[,p+1] + rnorm(n, 0.5))
  colnames(data) <- paste('X', 1:dim(data)[2], sep = "")
  colnames(data)[2*p+1] <- 'Y'
  return(data)
}

# 2.

rho = 0.1

# n >> p
n = 60
p = 5

donnees <- data.frame(CreateData(p, n, rho))
modelgd <- lm(Y ~ ., data = donnees)
summary(modelgd)

# Les variables sélectionnées sont la 1 et la 6, cad celles qui 
# ont été utilisées pour construire Y

# n > 2p
n = 15
p = 5

donnees <- data.frame(CreateData(p, n, rho))
modelgd <- lm(Y ~ ., data = donnees)
summary(modelgd)

# L'algorithme converge moins souvent qu'avant.

# n ≤ 2p
n = 5
p = 7

donnees <- data.frame(CreateData(p, n, rho))
modelgd <- lm(Y ~ ., data = donnees)
summary(modelgd)

# L'algorithme ne peut converger.

# 3.

rho = 0.9

# n >> p
n = 60
p = 5

donnees <- data.frame(CreateData(p, n, rho))
modelgd <- lm(Y ~ ., data = donnees)
summary(modelgd)

# Les variables sélectionnées sont la 1 et la 6, cad celles qui 
# ont été utilisées pour construire Y

# n > 2p
n = 15
p = 5

donnees <- data.frame(CreateData(p, n, rho))
modelgd <- lm(Y ~ ., data = donnees)
summary(modelgd)

# On retrouve parfois les bonnes variables après plusieurs essais, mais pas
# systématiquement, et on sélectionne parfois une des mauvaises variables en 
# raison de la forte colinéarité.

# n ≤ 2p
n = 5
p = 10

donnees <- data.frame(CreateData(p, n, rho))
modelgd <- lm(Y ~ ., data = donnees)
summary(modelgd)

# L'algorithme ne peut converger.

# Exercice 3 : Sélection de variables

# 1.

chen <- read.table("chenilles.txt",header=TRUE)
attach(chen)

# 2.

modelchen <- lm(NbNids ~ ., data = chen)
summary(modelchen)

# Les variables Altitude et Pente sont les plus importantes.

# 3.

reschen <- lm(NbNids ~ Altitude + Pente + NbPins + Hauteur + 
                       Diametre + Orient + NbStrat + Melange, data = chen)
summary(reschen)

# On constate que le R^2 (Adjusted R-squared) est plus grand pour le modèle 
# réduit que pour le modèle complet, il semble donc meilleur, ce que confirme 
# le test de Fisher suivant :

Ftest <- anova(reschen, modelchen)
Ftest

# Ici, la p-valeur vaut 0.98, le test est donc non significatif : 
# on ne rejette pas (H0), on garde le modèle réduit. En fait, pour trouver 
# le meilleur modèle, il faudrait faire p! tests (ce qu’on ne fera évidemment pas).

# 4.

summary(modelchen)

# on enlève du modèle la variable qui a la p-valeur la plus grande : Densite.

m10 <- lm(NbNids ~ . - Densite, data = chen)
summary(m10)

Ftest <- anova(modelchen, m10)
Ftest

# test non significatif, on préfère m10
# on enlève du modèle m10 la variable qui a la p-valeur la plus grande : HautMax.

m9 <- lm(NbNids ~ . -Densite -HautMax, data = chen)
summary(m9)

Ftest <- anova(m10, m9)
Ftest

# test non significatif, on préfère m9. On continue ainsi de suite jusqu'à m5.

m6 <- lm(NbNids ~ . -Densite -HautMax -Orient -Melange -NbPins, data = chen)
m5 <- lm(NbNids ~ . -Densite -HautMax -Orient -Melange -NbPins -NbStrat, data = chen)
summary(m5)

# tout est significatif, on s'arrête

Ftest <- anova(m6, m5)
Ftest

# Si on continue une fois de plus, le test deviendra significatif.
# Le modèle final retenu est le modèle m5 pour lequel il ne reste que 
# les variables Altitude, Pente, Hauteur et Diametre.

# 5.

step(modelchen, direction = "backward")

# on garde le modèle qui a la AIC-valeur la plus petite : -37.17.

m10 <- lm(NbNids ~ Altitude + Pente + Hauteur + Diametre + NbStrat, data = chen)
step(m10)

# Le modèle choisi est m6 (pour lequel la variable NbStrat a été ajoutée). 
# Ce n’est donc pas le même modèle que plus haut.