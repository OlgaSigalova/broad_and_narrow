library(tidyverse)
library(magrittr)

data(iris)


iris

names(iris)

iris$Species

iris[2,4]

iris[2, ]

iris %>%
  group_by(Species)  %>%
  summarize()


iris %>%
  summarize(mean(Sepal.Length),
            max(Sepal.Length),
            min(Sepal.Width))


iris %>%
  group_by(Species) %>%
  summarize(mean(Sepal.Length),
            max(Sepal.Length),
            min(Sepal.Width))


iris %>%
  filter(Sepal.Length>5)

# & - and,  | - or

iris %>%
  filter(Sepal.Length>5, Petal.Width>0.3)

iris %>%
  filter(grepl("set", Species))


iris %>%
  mutate(new_var = Sepal.Length + Sepal.Width)


ggplot(iris, aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point() +
  geom_smooth(method = "lm")


ggplot(iris, aes(x = Sepal.Length, fill = Species)) +
  geom_histogram(bins = 10, alpha = 0.5)


