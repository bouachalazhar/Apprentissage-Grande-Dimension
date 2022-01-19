# TD4 : Sélection de variables par pénalisation

# Exercice 1 : Découverte

## Importation des données


ozone <- read.table("/Users/bouacha_lazhar/OneDrive/Master MMA/M2 MMA/M2 S3/Apprentissage en Grande Dimension/Partie II/TP4  Sélection de variables par pénalisation-20201212/Ozone.txt")
dim(ozone)


### 1.


ozone2 <- ozone[,1:11]
model <- lm(maxO3v ~ ., data = ozone2)
summary(model)


### 2.


library(glmnet)

ozone2 <- scale(ozone2)
lambda_seq <- seq(0, 1, by = 0.001)
model_ridge <- glmnet(ozone2[,2:11], ozone2[,1], alpha = 0, lambda = lambda_seq, intercept = F)
plot(model_ridge)

library(ggplot2)

df = data.frame(lambda = rep(model_ridge$lambda, ncol(ozone2[,2:11])),
                beta = as.vector(t(model_ridge$beta)),
                variable = rep(colnames(ozone2[,2:11]), each = length(model_ridge$lambda)))
g1 = ggplot(df, aes(x = lambda, y = beta, col = variable)) + geom_line() +
  theme(legend.position = "bottom")
g1


### 3.


model_lasso <- glmnet(ozone2[,2:11], ozone2[,1], alpha = 1, lambda = lambda_seq, intercept = F)
plot(model_lasso)

df = data.frame(lambda = rep(model_lasso$lambda, ncol(ozone2[,2:11])),
                beta = as.vector(t(model_lasso$beta)),
                variable = rep(colnames(ozone2[,2:11]), each = length(model_lasso$lambda)))
g2 = ggplot(df, aes(x = lambda, y = beta, col = variable)) + geom_line() +
  theme(legend.position = "bottom")
g2


### 4.


model_EN <- glmnet(ozone2[,2:11], ozone2[,1], alpha = 0.5, lambda = lambda_seq, intercept = F)
plot(model_EN)

df = data.frame(lambda = rep(model_EN$lambda, ncol(ozone2[,2:11])),
                beta = as.vector(t(model_EN$beta)),
                variable = rep(colnames(ozone2[,2:11]), each = length(model_EN$lambda)))
g3 = ggplot(df, aes(x = lambda, y = beta, col = variable)) + geom_line() +
  theme(legend.position = "bottom")
g3


### 5.

# ridge
ridge_cv <- cv.glmnet(ozone2[,2:11], ozone2[,1], alpha = 0, lambda = lambda_seq,
                      nfolds = 10, intercept = F)
best_lambda <- ridge_cv$lambda.min
plot(ridge_cv)
ridge <- glmnet(ozone2[,2:11], ozone2[,1], alpha = 0, lambda = best_lambda, intercept = F)
ridge$beta

# lasso
lasso_cv <- cv.glmnet(ozone2[,2:11], ozone2[,1], alpha = 1, lambda = lambda_seq, nfolds = 10,
                      intercept = F)
best_lambda <- lasso_cv$lambda.min
plot(lasso_cv)
lasso <- glmnet(ozone2[,2:11], ozone2[,1], alpha = 1, lambda = best_lambda, intercept = F)
lasso$beta

# elastic net
en_cv <- cv.glmnet(ozone2[,2:11], ozone2[,1], alpha = 0.5, lambda = lambda_seq, nfolds = 10,
                   intercept = F)
best_lambda <- en_cv$lambda.min
plot(en_cv)
en <- glmnet(ozone2[,2:11], ozone2[,1], alpha = 0.5, lambda = best_lambda, intercept = F)
en$beta

# graphique qui résume
regusuel <- lm(ozone2[,1]~ozone2[,2:11]-1) # modèle linéaire sans intercept
df4 <- data.frame(x=rep(colnames(ozone2[,2:11]),4),
                  coef=c(as.vector(regusuel$coefficients),
                         as.vector(coef(ridge_cv,s=ridge_cv$lambda.min)[-1]),
                         as.vector(coef(lasso_cv)[-1]),as.vector(coef(en_cv)[-1])),
                  reg=c(rep("reg.lin",ncol(ozone2[,2:11])),
                        rep("ridge",ncol(ozone2[,2:11])),
                        rep("lasso",ncol(ozone2[,2:11])),
                        rep("ElasticNet",ncol(ozone2[,2:11]))))
g4 <- ggplot(df4)+ geom_point(aes(x=x,y=coef,col=reg))
g4


# Exercice 2 : Comparaison des méthodes

## Jeu de données 1 : petit signal et beaucoup de bruit


p <- 5000
n <- 1000
real_p <- 15
x <- matrix(rnorm(n*p), nrow=n, ncol=p)
y <- apply(x[,1:real_p], 1, sum) + rnorm(n)

train_rows <- sample(1:n, .66*n)
x.train1 <- x[train_rows,]
x.test1 <- x[-train_rows,]
y.train1 <- y[train_rows]
y.test1 <- y[-train_rows]


## Jeu de données 2 : gros signal et beaucoup de bruit


p <- 5000
n <- 1000
real_p <- 1000
x <- matrix(rnorm(n*p), nrow=n, ncol=p)
y <- apply(x[,1:real_p], 1, sum) + rnorm(n)

train_rows <- sample(1:n, .66*n)
x.train2 <- x[train_rows,]
x.test2 <- x[-train_rows,]
y.train2 <- y[train_rows]
y.test2 <- y[-train_rows]


## Jeu de données 3 : signal varié et variables corrélées

library(MASS)
p <- 50
n <- 100
CovMatrix <- outer(1:p, 1:p, function(x,y) {.7^abs(x-y)})
x <- mvrnorm(n, rep(0,p), CovMatrix)
y <- 10 * apply(x[, 1:2], 1, sum) + 5 * apply(x[, 3:4], 1, sum) + 
  apply(x[, 4:14], 1, sum) + rnorm(n)

train_rows <- sample(1:n, .66*n)
x.train3 <- x[train_rows,]
x.test3 <- x[-train_rows,]
y.train3 <- y[train_rows]
y.test3 <- y[-train_rows]


### 1.

# Jeu 1
fit.lasso <- glmnet(x.train1, y.train1, family="gaussian", alpha=1)
fit.ridge <- glmnet(x.train1, y.train1, family="gaussian", alpha=0)
fit.elnet <- glmnet(x.train1, y.train1, family="gaussian", alpha=.5)

# Jeu 2
fit.lasso <- glmnet(x.train2, y.train2, family="gaussian", alpha=1)
fit.ridge <- glmnet(x.train2, y.train2, family="gaussian", alpha=0)
fit.elnet <- glmnet(x.train2, y.train2, family="gaussian", alpha=.5)

# Jeu 3
fit.lasso <- glmnet(x.train3, y.train3, family="gaussian", alpha=1)
fit.ridge <- glmnet(x.train3, y.train3, family="gaussian", alpha=0)
fit.elnet <- glmnet(x.train3, y.train3, family="gaussian", alpha=.5)

### 2.

plot(fit.lasso)
plot(fit.ridge)
plot(fit.elnet)

### 3.


# Jeu 1
lasso <- cv.glmnet(x.train1, y.train1, type.measure="mse",alpha=1,family="gaussian")
EN <- cv.glmnet(x.train1, y.train1, type.measure="mse",alpha=0.5,family="gaussian")
ridge <- cv.glmnet(x.train1, y.train1, type.measure="mse",alpha=0,family="gaussian")

lasso.model <- glmnet(x.train1, y.train1, alpha=1, family="gaussian", lambda=lasso$lambda.min)
EN.model <- glmnet(x.train1, y.train1, alpha=0.5, family="gaussian", lambda=EN$lambda.min)
ridge.model <- glmnet(x.train1, y.train1, alpha=0, family="gaussian", lambda=ridge$lambda.min)
yhat_lasso <- predict(lasso.model, s=lasso.model$lambda.min, newx=x.test1)
yhat_EN <- predict(EN.model, s=lasso.model$lambda.min, newx=x.test1)
yhat_ridge <- predict(ridge.model, s=lasso.model$lambda.min, newx=x.test1)
mse_lasso <- mean((y.test1 - yhat_lasso)^2)
mse_EN <- mean((y.test1 - yhat_EN)^2)
mse_ridge <- mean((y.test1 - yhat_ridge)^2)
c(mse_lasso,mse_EN,mse_ridge) 


# lasso est le meilleur (normal, il est fait pour ça)

### 4.

# jeu 2: ridge est le meilleur (bon en prédiction mais pas vraiment interprétable)
# jeu 3: EN est le meilleur


# Exercice 3 : Comparaison pls et Lasso

## Importation des données


library(plsgenomics)

data(Colon)
length(Colon)
X <- Colon$X
Y <- Colon$Y
Y <- Y-1
gene <- Colon$gene.names
Colon <- data.frame(X=I(X),Y=Y)

### 1.


train <- rbinom(length(Colon$Y),1,2/3)
Colon.train <- c()
Colon.test <- c()
Colon.train$X <- Colon$X[train==1,]
Colon.train$Y <- Colon$Y[train==1]
Colon.test$X <- Colon$X[train==0,]
Colon.test$Y <- Colon$Y[train==0]


### 2.


colon.glm <- glmnet(Colon.train$X, Colon.train$Y, family="binomial", alpha = 0.5)
plot(colon.glm) # elastic net
colon.glm <- glmnet(Colon.train$X,Colon.train$Y,family="binomial",alpha=1)
plot(colon.glm) # ici, il s'agit du lasso


### 3.

colon.cv <- cv.glmnet(Colon.train$X,Colon.train$Y,nfolds=5)
colon.glmnet.model <- glmnet(Colon.train$X,Colon.train$Y,family="binomial",nlambda=1,lambda=colon.cv$lambda.min)
length(which(abs(colon.glmnet.model$beta)>0))


### 4.

library(pls)

plscolon <- plsr(Colon.train$Y~Colon.train$X,ncomp=30,scale=TRUE,validation="CV",segments=5)
summary(plscolon)
ncomponents=3 # 75% de la variance des Y expliquée
pls.reduction.matrix <- plscolon$loadings[,1:ncomponents]
Colon.train$pls.reducedX <- Colon.train$X %*% pls.reduction.matrix
Colon.test$pls.reducedX <- Colon.test$X %*% pls.reduction.matrix
pls.model <- glm(Y~pls.reducedX,data=Colon.train,family=binomial(link="logit"))
test.prediction <- predict(pls.model,newdata = Colon.test,type="response")
rbind(test.prediction,Colon.test$Y) # 4 erreurs uniquement


### 5.


colonpredict.glmnet <- predict.glmnet(colon.glmnet.model,Colon.test$X,type="response")
sum((1/length(Colon.test$Y))*(colonpredict.glmnet-Colon.test$Y)^2)
sum((1/length(Colon.test$Y))*(test.prediction-Colon.test$Y)^2)