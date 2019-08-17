
options(scipen = 999)

p <- ggplot(df1, aes(x = interquantile_width)) +
  geom_density()


plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

set.seed(1)
wait <- df1$interquantile_width
mixmdl <- normalmixEM(wait, k = 2)


data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.5) +
  ylab("count")

post.df <- as.data.frame(cbind(x = mixmdl$x, mixmdl$posterior)) %>%
  unique()
head(post.df, 10)