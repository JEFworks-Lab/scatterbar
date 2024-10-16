# Define the cities
cities <- c("New York", "Los Angeles", "Mexico City", "Buenos Aires",
            "London", "Cairo", "Beijing", "Sydney", "Tokyo")

# Function to generate random proportions that add up to 1
generate_proportions <- function() {
  props <- runif(4)
  props / sum(props)
}

# Create the data frame with random proportions
city_data <- t(sapply(1:length(cities), function(x) generate_proportions()))
colnames(city_data) <- c("Dogs", "Cats", "Fish", "Birds")
df <- data.frame(city_data)
rownames(df)<- cities
pos_cities <- data.frame(Cities = cities,
  x = c(8, 1, 4, 7, 11, 11, 16, 20, 20),
  y = c(19, 17, 14, 2, 18, 11, 9, 4, 20)
)
row.names(pos_cities) <- pos_cities$Cities
pos_cities <- pos_cities %>% select(-Cities)

scatterbar::scatterbar(df, pos_cities, plot_title="Pet Owners Across the World")
