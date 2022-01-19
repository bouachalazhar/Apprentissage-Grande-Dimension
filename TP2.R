# TD2 : Sélection par tests multiples

# Exercice 1 : Tests Multiples sur un jeu de données simulées

## Description du jeu de données


load("/Users/bouacha_lazhar/OneDrive/Master MMA/M2 MMA/M2 S3/Apprentissage en Grande Dimension/Partie II/TP2  Tests multiples-20201116/genes.Rdata")
data <- genes$data
condition <- genes$condition
statut <- genes$statut


### 1.

# On peut faire un test de comparaison de moyenne avec échantillons appariés :


t.test(data[condition == 1, 1], data[condition == 2, 1])


# La p-valeur > 5%, le test est donc non significatif. Le gène 1 est non différentiellement exprimé.

### 2.


pval <- apply(data, 2, function(x) t.test(x[condition == 1], x[condition == 2])$p.value)
length(which(pval < 0.05)) # ceux qui sont détéctés positifs


# Sur les 5000 gènes, 334 sont identifiés comme différentiellement exprimés. En tout (voir code suivant), 96 gènes sont réellement différentiellement exprimés.


length(which(statut != 0)) # ceux qui sont réellement positifs


# Pour connaître le nombre de faux positifs, on compare la liste de gènes identifiés comme différentiellement exprimés avec ceux qui le sont réellement :


FP <- length(intersect(which(statut == 0), which(pval < 0.05)))
FP


# Le nombre de Faux Positifs est très élevé mais cela n’est pas surprenant puisque environ 5% des gènes non différentiellement exprimés sont classés à tort. Cela représente environ (5000 - 100)*0.05 = 245 gènes, où 100 : Vrais Positifs

### 3.


length(which(pval < 0.05/5000))
statut[which(pval < 0.05/5000)]
FP <- length(intersect(which(statut == 0), which(pval < 0.05/5000)))
FP

# Si on applique la procédure de Bonferroni, on ne garde que 20 gènes (procédure très conservative). Il n’y a aucun faux positif.

### 4.


length(which(pval < (1-(1-0.05)^(1/5000))))
statut[which(pval < (1-(1-0.05)^(1/5000)))]
FP <- length(intersect(which(statut == 0), which(pval < (1-(1-0.05)^(1/5000)))))
FP


# Si on applique la procédure de Sidàk, on ne garde que 20 gènes identique à Bonferroni. Il n’y a aucun faux positif.

### 5.


holmbonferroni <- function(pval, alpha){ # création d'une fonction holm-bonferroni
  m <- length(pval)
  # 1. ordonne les p-valeurs
  sortedpval <- sort(pval)
  # 2. détermine
  x <- which(sortedpval <= alpha/(m+1-c(1:m)))
  I <- max(x) # rang à partir duquel il faut couper la liste qu'on renvoie
  # 3. extraire
  selected = c()
  if (I > 1){
    selected = names(sortedpval)[1:I]
  }
  return(selected)
}

holmbonferroni(pval, 0.05)
length(holmbonferroni(pval, 0.05))


# Les résultats sont une nouvelle fois les mêmes.

### 6.


benjaminihochberg <- function(pval, alpha){ # création d'une fonction Benjamini-Hochberg
  m <- length(pval)
  # 1. ordonne les p-valeurs
  sortedpval <- sort(pval)
  # 2. détermine
  x <- which(sortedpval <= alpha*(c(1:m)/m))
  I <- max(x) # rang à partir duquel il faut couper la liste qu'on renvoie
  # 3. extraire
  selected = c()
  if (I > 0){
    selected = names(sortedpval)[1:I]
  }
  return(selected)
}

benjaminihochberg(pval, 0.05)
length(benjaminihochberg(pval, 0.05))
FP <- length(intersect(names(which(statut == 0)), benjaminihochberg(pval,.05)))
Positifs <- length(benjaminihochberg(pval,.05))
FDP <- FP/Positifs
FDP


# La procédure de Benjamini-Hochberg est moins conservative. On sélectionne plus de gènes, quitte à commettre un peu plus d’erreurs (3 faux positifs ici). Le taux de faux positifs (FDP) est de l’ordre de 5%.

### 7.


padj <- p.adjust(pval, method = "BH")
which(padj < 0.05) # coincide bien
length(which(padj < 0.05))


# On peut essayer avec hochberg, bonferroni


adj <- p.adjust(pval, method = "hochberg")
which(padj < 0.05) # coincide bien
length(which(padj < 0.05))

adj <- p.adjust(pval, method = "bonferroni")
which(padj < 0.05) # coincide bien
length(which(padj < 0.05))


# Exercice 2 : Tests Multiples sur un jeu de données réelles

## Chargement et description du jeu de données 

#BiocManager::install("golubEsets")
library(golubEsets)
data("Golub_Merge")
load("/Users/bouacha_lazhar/OneDrive/Master MMA/M2 MMA/M2 S3/Apprentissage en Grande Dimension/Partie II/TP2  Tests multiples-20201116/Exp.Rdata")
Exp <- exprs(Golub_Merge)
dim(Exp)

### 1.


hist(Exp, main="")

Exp[Exp < 100] <- 100
Exp[Exp > 16000] <- 16000

emax <- apply(Exp, 1, max)
emin <- apply(Exp, 1, min)
R <- which(emax/emin > 5)
S <- which(emax - emin > 500)
Exp <- Exp[intersect(R, S),]
dim(Exp)

Exp <- log(Exp,10)
Exp <- scale(t(Exp), center=TRUE, scale=TRUE) # ne pas oublier de transposer la matrice

hist(Exp, main="")


### 2.


load("/Users/bouacha_lazhar/OneDrive/Master MMA/M2 MMA/M2 S3/Apprentissage en Grande Dimension/Partie II/TP2  Tests multiples-20201116/Clin.Rdata")
ALL <- which(Golub_Merge$ALL.AML == "ALL")
AML <- which(Golub_Merge$ALL.AML == "AML")

pval <- apply(Exp, 2, function(x){
  t.test(x[ALL], x[AML])$p.value
})

pvalBH <- p.adjust(pval,method="BH")
length(which(pvalBH<0.01))


# Sur les 3571 gènes, 571 sont identifiés comme différentiellement exprimés avec un contrôle de la FDR de 1%.


DEgenes <- Exp[,which(pvalBH<0.01)]
dim(DEgenes)


### 3.


library(gplots)

rownames(DEgenes) <- Golub_Merge$ALL.AML
heatmap(DEgenes, col=greenred(75), ylab="Conditions", xlab="Genes")


# Des groupes de gènes apparaissent distinctement et pourraient permettre de différencier les deux types de leucémie.

### 4.


bestgenes <- order(pvalBH)[1:50]
DEbest <- Exp[,bestgenes]
dim(DEbest)

rownames(DEbest) <- Golub_Merge$ALL.AML
heatmap(DEbest, col=greenred(75), ylab="Conditions", xlab="Genes")

library(randomForest)

I <- rbinom(nrow(DEbest),1,.8) # pour créer les training et test sets
Train <- which(I==1)
Test <- which(I==0)
rf <- randomForest(x = DEbest[Train,],y=Golub_Merge$ALL.AML[Train],
                   xtest=DEbest[Test,], ytest = Golub_Merge$ALL.AML[Test])
rf

# La méthode de classification fonctionne très bien, aucune erreur n’est commise.