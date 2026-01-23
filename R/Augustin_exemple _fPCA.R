library(fda)
library(tidyverse)

f <- funfdaf <- function(x){
  r1 <- rnorm(1)
  r2 <- rnorm(1)
  return(x^(2+r1)+3*r2*x)
}


df <- data.frame(ind = paste0("I",0:9),
                 t(as.matrix(sapply(X = 1:10, FUN = function(x){as.matrix(sapply(X = 1:5, FUN = f))})))
                 )
colnames(df) <- c("ind",1:5)

df %>% pivot_longer(-ind,
                    names_to = "x",
                    values_to = "y") %>% 
  group_by(ind) %>% 
  ggplot()+ 
  aes(x = x,
      y = y, col = ind) +
  geom_line(aes(group = ind))



x = 1:5
xrange = c(1,5)
nbsplines = 5

basis <- create.bspline.basis(rangeval = xrange, nbasis = nbsplines)
# plot(basis)
fdParobj <- fdPar(basis, Lfdobj = 2, lambda = 1e-4)
fdobj <- smooth.basis(argvals = x, y = t(as.matrix(df[-1])), fdParobj = fdParobj)$fd


# df_smooth <- eval.fd(evalarg = seq(from = 1, to = 5, by = 0.01), fdobj = fdobj)
# 
# 
# as.data.frame(df_smooth) %>% mutate(x = seq(from = 1, to = 5, by = 0.01)) %>%
#   pivot_longer(cols = -x, names_to = "ind", values_to = "y") %>%
#   group_by(ind) %>%
#   ggplot(aes(x = x, y = y, col = ind)) +
#   geom_line(aes(group = ind))


res.fpca <- pca.fd(fdobj, nharm = 5)

coord_ind <- res.fpca$scores
# plot(coord_ind[,1],coord_ind[,2])
beta0 <- eval.fd(evalarg = seq(from = 1, to = 5, by = 0.01), res.fpca$meanfd)

harm <- eval.fd(evalarg = seq(from = 1, to = 5, by = 0.01), res.fpca$harmonics)

courbes_1 = beta0 

for(i in 1:dim(coord_ind)[2]){
  courbes_1 <- courbes_1 + coord_ind[1,i]*harm[,i]
}

# Affichage des individus
data.frame(coord_ind[,1:2]) %>% 
  ggplot() +
  aes(x = X1, y = X2) +
  geom_point() +
  geom_hline(yintercept = 0, col = "grey30", linetype = "dashed") +
  geom_vline(xintercept = 0, col = "grey30", linetype = "dashed") +
  theme_minimal() +
  ggsave("R/Images/Nuages_ind_exemple.svg")


# Recomposistion des courbes des individus
data.frame(x = seq(from = 1, to = 5, by = 0.01), c_reel = df_smooth[,1], c_mod = courbes_1) %>% 
  pivot_longer(-x,
               values_to = "y",
               names_to = "obj") %>% 
  group_by(obj) %>% 
  ggplot() +
  aes(x = x, y = y, col = obj) +
  geom_line(aes(group = obj))

# Affichage fonction moyenne
data.frame(x = seq(from = 1, to = 5, by = 0.01), mean = beta0) %>% 
  ggplot() +
  aes(x = x, y = mean) +
  geom_line() +
  theme_minimal() +
  ggsave("R/Images/fonction_mean_exemple.svg")


# Affichage des harmoniques
data.frame(x = seq(from = 1, to = 5, by = 0.01), harm[,1:3]) %>% 
  pivot_longer(-x,
               names_to = "group",
               values_to = "y") %>% 
  group_by(group) %>% 
  ggplot() +
  aes(x = x, y = y, col = group) +
  geom_line(aes(group = group)) +
  theme_minimal() +
  ggsave("R/Images/fonction_harm_exemple.svg")


cumsum(res.fpca$varprop)







