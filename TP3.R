# TD3 : Réduction de dimension par ACP et PLS

# Exercice 1 : Interprétation des graphiques de l’ACP

## Importation des données
options(ggrepel.max.overlaps = Inf)
body <- read.table("/Users/bouacha_lazhar/OneDrive/Master MMA/M2 MMA/M2 S3/Apprentissage en Grande Dimension/Partie II/TP3  Réduction de dimension par ACP et PLS-20201212/body_full.csv"
                   ,header=TRUE, sep=";",dec=",")
dim(body)
library(FactoMineR)

### 1.


boxplot(body[,-25],las=3,cex.axis=0.65)


# Ici, on voit que les données sont très disparates, il vaut sans doute mieux les normaliser.


body.norm <- scale(body[,-25])
boxplot(body.norm,las=3,cex.axis=0.65)


### 2.


res.PCA.body.full <- PCA(body.norm)


# On peut vouloir garder l’ensemble des variables (y compris qualitatives) pour effectuer l’ACP. Auquel cas, on utilise la fonction PCA() sur le jeu de données global en précisant les variables qualitatives (ici, la variable 25).


res.PCA.body.full <- PCA(body,quali.sup=25)


### 3.


plot(res.PCA.body.full, choix="ind", habillage=25, label="none")


# Les deux composantes principales permettent de différencier les hommes des femmes.

### 4.


plot(res.PCA.body.full, choix="var")


# Le graphique est relativement complexe pour trouver les variables qui définissent les composantes principales.
# On voit cependant que :
#  • le 1er axe est défini par la largeur ou carrure des individus (waist.girth, weight, ankle.girth, chest.depth, chest.girth), correspondant ainsi plutôt aux hommes,
#  • le 2ième axe est défini par l’épaisseur des individus (thigh girth, hip.length), correspondant ainsi plutôt aux femmes.


# Exercice 2 : Problème de normalisation

## Chargement des données


athlete <- read.csv("/Users/bouacha_lazhar/OneDrive/Master MMA/M2 MMA/M2 S3/Apprentissage en Grande Dimension/Partie II/TP3  Réduction de dimension par ACP et PLS-20201212/athle_records.csv"
                    ,sep = "\t",header=TRUE,dec=",")
rownames(athlete) <- athlete$X
athlete <- athlete[,-1]
rownames(athlete)[3] <- "Bresil"
rownames(athlete)[14] <- "Jamaique"
rownames(athlete)[18] <- "NouvelleZelande"
rownames(athlete)[23] <- "Suede"
dim(athlete)


### 1.


boxplot(athlete,las=2,cex.axis=0.7)


# La variable marathon risque d’avoir trop de poids.

### 2.


res <- prcomp(athlete)
res


# La composante 1 semble être uniquement déterminée par marathon et semi-marathon.

### 3.


biplot(res)


# On a effectivement l’impression que la composante principale 1 est entièrement définie par la variable marathon. Il n’y a pas grand chose à dire concernant la 2ième composante principale.

### 4.


plot(res)


# On se rend compte que la composante 1, c’est-à-dire la variable marathon suffit à expliquer la quasi-totalité des données.

### 5.


athlete.log <- log(athlete)
boxplot(athlete.log)
res.log <- prcomp(athlete.log)
plot(res.log)


# L’éboulis des valeurs propres et la règle du coude nous indique qu’il faudrait garder deux composantes principales.


biplot(res.log)


# L’ACP sur le jeu de données normalisé est beaucoup plus pertinent. Ici, on voit que les composantes principales séparent les pays performants en épreuve d’endurance (marathon, semi-marathon - composante 1) de ceux performants en distance courte (100m, 200m - composante 2). On retrouve la Jamaïque sur la composante 2, l’Ethiopie et le Kenya sur la composante 1 (faits connus).

### 6.


biplot(res.log,choices=2:3)


# La composante 3 sépare encore les pays selon d’autres types d’épreuves (marathon, semi-marathon vs 800m, 1500m). Attention cependant, la part de variance expliquée par la composante 3 est relativement faible.

### 7.


PCA(athlete)


# Les résultats sont identiques.

# Exercice 3 : Application à un jeu de données biologiques

## Importation des données


load("/Users/bouacha_lazhar/OneDrive/Master MMA/M2 MMA/M2 S3/Apprentissage en Grande Dimension/Partie II/TP3  Réduction de dimension par ACP et PLS-20201212/Prostate.Rdata")


### 1.


boxplot(Prostate, las=1,cex.axis=0.75)
Pros <- scale(Prostate)
boxplot(Pros, las=1,cex.axis=0.75)


### 2.


res.pca <- PCA(Pros)


### 3.


eigenvalues <- res.pca$eig
eigenvalues

barplot(eigenvalues[, 2], names.arg=1:nrow(eigenvalues),
        main = "Eboulis des valeurs propres",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")

lines(x = 1:nrow(eigenvalues), eigenvalues[, 2],
      type="b", pch=19, col = "red")


# Suivant le critère du coude, on ne garderait que les deux 1ers axes. Ici, le 1er axe retient 43% de l’inertie totale, le 2ième 18% (pour un total de 60%).

### 4.


PCA(Pros)


# Les graphiques sont assez clairs : le 1er axe est défini par les variables lpsa, lcavol, gleason, pgg45, svi et lcp, le 2ième par les variables lbph, lweight et age. Un peu de recherche bibliographique nous aiderait certainement à mieux comprendre à quoi correspondent ces axes.

### 5.

library(factoextra)

fviz_contrib(res.pca, choice = "ind", axes = 1,top=40)
fviz_contrib(res.pca, choice = "ind", axes = 2,top=40)
fviz_contrib(res.pca, choice = "var", axes = 1)
fviz_contrib(res.pca, choice = "var", axes = 2)


# Ces graphiques offrent une visualisation différente de la contribution des variables-individus par composante.

# Exercice 4 : ACP-PLS sur un jeu de données génomique

## Chargement des données


library(plsgenomics)

data(Colon)
length(Colon)
X <- Colon$X
Y <- Colon$Y
Y <- Y-1
gene <- Colon$gene.names
colon <- data.frame(X=I(X),Y=Y)

### 1.


train <- rbinom(length(colon$Y),1,2/3)
colon.train <- c()
colon.test <- c()
colon.train$X <- colon$X[train==1,]
colon.train$Y <- colon$Y[train==1]
colon.test$X <- colon$X[train==0,]
colon.test$Y <- colon$Y[train==0]


# Soit Y la variable binaire indiquant le tissu d’origine des échantillons et X le jeu de données de taille n × p
# contenant les données d’expression des p = 2000 gènes sur les n échantillons. On met alors en place le modèle
# de régression logistique suivant :
#  log(P(Y = 1)/P(Y = 0)) = X*beta,
# où beta (à estimer) indique le lien existant entre Y et X. Sur R, cela donne :


model <- glm(Y~X, family="binomial", data=colon.train)
summary(model)


# Le modèle de régression logistique ne fonctionne pas puisque nous sommes dans un cadre de grande dimension
# p > n.

### 2.


library(pls)

pcrcolon <- pcr(Y ~ X, data = colon.train, ncomp=30,scale=FALSE,validation="none")
summary(pcrcolon)


# Il faudrait garder 6 composantes pour expliquer 75% de la variance du nuage.

### 3.


ncomponents <- 6
reduction.matrix <- pcrcolon$loadings[,1:ncomponents] # matrice de l'ACP réduite aux 6 premières composantes
colon.train$reducedX <- colon.train$X %*% reduction.matrix
colon.test$reducedX <- colon.test$X %*% reduction.matrix


### 4.


pca.model <- glm(Y~reducedX,data=colon.train,family=binomial(link="logit"))
test.prediction <- predict(pca.model,newdata = colon.test,type="response")
test.prediction.bin <- test.prediction
test.prediction.bin[test.prediction.bin > 0.5] <- 1
test.prediction.bin[test.prediction.bin <= 0.5] <- 0
rbind(test.prediction.bin,colon.test$Y) # 4 mal classés sur 19
sum((test.prediction-colon.test$Y)^2)


### 5.


plscolon <- plsr(Y~X,data=colon.train,ncomp=30,scale=TRUE,validation="CV",segments=5)
summary(plscolon)

ncomponents=4 # 75% de la variance des Y expliquée
pls.reduction.matrix <- plscolon$loadings[,1:ncomponents]
colon.train$pls.reducedX <- colon.train$X %*% pls.reduction.matrix
colon.test$pls.reducedX <- colon.test$X %*% pls.reduction.matrix

pls.model <- glm(Y~pls.reducedX,data=colon.train,family=binomial(link="logit"))
test.prediction <- predict(pls.model,newdata = colon.test,type="response")
rbind(test.prediction,colon.test$Y) # 2 erreurs uniquement
sum((test.prediction-colon.test$Y)^2)


# Inutile de binariser manuellement la prédiction.