
fixed <- fixed
moving <- moving

Na <- 1000
Nb <- 1000

# Sample from fixed and moving images for density estimation

coords <- sample_coords(dim(fixed), size=Na)
sampA_f <- intensities_from_coords(fixed, coords, plot=FALSE, new=TRUE)

sampA_m <- intensities_from_coords(moving, coords, plot=FALSE, new=TRUE)




