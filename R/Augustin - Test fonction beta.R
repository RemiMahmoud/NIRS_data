

pheno_var <- "TSW"

# Extraire la matrice NIRS (lignes = individus, colonnes = longueurs d'onde)
X_mat_train <- as.matrix(full_data_train[, grep("^[0-9]+$", colnames(full_data_train))])
X_mat_test <- as.matrix(full_data_test[, grep("^[0-9]+$", colnames(full_data_test))])
wn <- as.numeric(colnames(X_mat_train))  # 817, ..., 2764

Y_train <- full_data_train[[pheno_var]]
Y_test <- full_data_test[[pheno_var]]

nbasis <- 25
basis <- create.bspline.basis(rangeval = range(wn), nbasis = nbasis)
fdPar_obj <- fdPar(basis, Lfdobj = 2, lambda = 1e-4)

fdobj <- smooth.basis(wn, t(X_mat_train), fdPar_obj)$fd

# Régression scalar-on-function
fReg <- fRegress(Y_train ~ nirs_fdobj, data = list(Y_train = Y_train, nirs_fdobj = fdobj))


fdobj_test = smooth.basis(wn, t(X_mat_test), fdPar_obj)$fd


# Coefficient fonctionnel
beta_fd <- fReg$betaestlist$nirs_fdobj$fd

# Intercept
beta0 <-  as.numeric(eval.fd(0,  fReg$betaestlist$const$fd))

# Y predit
yhat_pred_test <- beta0 + inprod(fdobj_test, beta_fd)[,1]

# RMSE
RMSE_test <- sqrt(mean((Y_test - yhat_pred_test)^2))

# Affichage fonction Beta
beta <- eval.fd(evalarg = wn, fdobj = beta_fd)

ggplot(data.frame(wn = wn, beta = beta)) +
  aes(x = wn, y = beta) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey20") +
  geom_vline(xintercept = wn[which.max(beta)], linetype = "dashed", col = "red", linewidth = 1) +
  theme_minimal() + 
  ggtitle("Fonction Beta d'une régression fonctionnelle sur TSW") +
  xlab("longueur d'ondes (nm)") +
  ylab("Beta") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"))

# Création des courbes
    # Intérêt de faire ça pour plot ? pas sûr... ==> Oui car il y a un lissage sur les courbes
fdobj_full_dataset <- smooth.basis(wn, t(as.matrix(full_data[,-c(1,1950:1961)])), fdPar_obj)$fd
X <- eval.fd(evalarg = wn, fdobj = fdobj_full_dataset)


# Affichage des courbes et de leur valeur de TSW

# On change les noms de colonnes pour avoire les valeurs de TSW et individu

new_names <- paste0("I",1:208,"_",pheno_VG_pred$TSW)
# new_names
colnames(X) <- new_names


# dessin des courbes
  # TSW continu
as.data.frame(X) %>% mutate(wn = wn) %>% pivot_longer(-c(wn),
                                                      names_to = "Individu",
                                                      values_to = "Abs") %>% 
  mutate(TSW = as.numeric(sapply(strsplit(Individu,"_"), `[`, 2)), Ind_Id = sapply(strsplit(Individu,"_"), `[`, 1)) %>%
  group_by(Ind_Id) %>% 
  ggplot() +
  aes(x = wn, y = Abs, col = TSW) +
  geom_line(aes(group = Individu)) +
  geom_vline(xintercept = wn[which.max(beta)], col = "red", linetype = "dashed") +
  scale_color_continuous(palette = c("red","white","blue")) +
  theme_minimal() + 
  ggtitle("Spectre d'absorbtion sur graines") +
  xlab("longueur d'ondes (nm)") +
  ylab("Absorbance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) #+
  # ggsave("R/Images/Beta.png", width = 11, height = 6)
  
  # Zoom sur le beta max
  # ylim(-0.005,0) +
  # xlim(1050,1150) +
  # ggsave("R/Images/Zoom_beta_max.png", width = 11, height = 8)
  
  
  # Zoom sur un beta nul
  # ylim(0,0.02) +
  # xlim(2650,2750) +
  # ggsave("R/Images/Zoom_beta_nul.png", width = 11, height = 8)
  
  

  # TSW discret
p1 <- as.data.frame(X) %>% mutate(wn = wn) %>% pivot_longer(-c(wn),
                                                   names_to = "Individu",
                                                   values_to = "Abs") %>% 
  mutate(TSW = as.numeric(sapply(strsplit(Individu,"_"), `[`, 2)), Ind_Id = sapply(strsplit(Individu,"_"), `[`, 1)) %>% 
  mutate(TSW_class = as.factor(floor(TSW))) %>% 
  group_by(Ind_Id) %>% 
  ggplot() +
  aes(x = wn, y = Abs, col = TSW_class) +
  geom_line(aes(group = Individu)) +
  geom_vline(xintercept = wn[which.max(beta)], col = "red", linetype = "dashed")+
  geom_vline(xintercept = wn[which.min(beta)], col = "blue", linetype = "dashed")


# Recherche des points "exceptionels"
library(plotly)
ggplotly(p1)

# Affichage des courbes beta, X1 et beta*X1

p1 <- data.frame(wn = wn, beta = beta, X = X[,1]) %>% mutate(prod = beta*X) %>% ggplot() +
  aes(x = wn, y = X) + 
  geom_line() +
  geom_vline(xintercept = wn[which.max(beta)], col = "red", linetype = "dashed") +
  geom_vline(xintercept = wn[which.min(beta)], col = "blue", linetype = "dashed") +
  theme_minimal() + 
  ggtitle("Spectre d'absorbance") +
  xlab("longueur d'ondes (nm)") +
  ylab("Absorbance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14))



p2 <- data.frame(wn = wn, beta = beta, X = X[,1]) %>% mutate(prod = beta*X) %>% ggplot() +
  aes(x = wn, y = beta) + 
  geom_line() +
  geom_vline(xintercept = wn[which.max(beta)], col = "red", linetype = "dashed")+
  geom_vline(xintercept = wn[which.min(beta)], col = "blue", linetype = "dashed") +
  theme_minimal() + 
  ggtitle("Fonction Beta (rf sur TSW)") +
  xlab("longueur d'ondes (nm)") +
  ylab("Beta") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14))

p3 <- data.frame(wn = wn, beta = beta, X = X[,1]) %>% mutate(prod = beta*X) %>% ggplot() +
  aes(x = wn, y = prod) + 
  geom_line() +
  geom_vline(xintercept = wn[which.max(beta)], col = "red", linetype = "dashed")+
  geom_vline(xintercept = wn[which.min(beta)], col = "blue", linetype = "dashed") +
  theme_minimal() + 
  ggtitle("Spectre * Beta") +
  xlab("longueur d'ondes (nm)") +
  ylab("Beta") +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14))

p1 + p2 + p3


# Evolution de l'intégrale (mais on s'en fiche)

list_inn_prod <- rep(0, times = length(wn))
list_inn_prod[1] <- beta[1] * X[1,1]
for(i in 2:length(wn)){
  list_inn_prod[i] <- list_inn_prod[i-1] + (beta[i] * X[i,1])
}


data.frame(wn = wn, inner_product = list_inn_prod) %>% ggplot() +
  aes(x = wn, y = inner_product) +
  geom_line() +
  geom_hline(yintercept = inprod(beta_fd, fdobj)[1], col = "red", linetype = "dashed") +
  geom_vline(xintercept = 1090, col = "red", linetype = "dashed")

