# RStan needs these headers when compiling a model at runtime, even when RStan
# itself came from a binary package. pak need not install LinkingTo-only packages
# for binaries, so the native CI job requests BH and RcppEigen explicitly.
headers <- c(
  BH = system.file("include", "boost", "version.hpp", package = "BH"),
  RcppEigen = system.file("include", "Eigen", "Core", package = "RcppEigen")
)
missing <- names(headers)[!nzchar(headers) | !file.exists(headers)]
if (length(missing)) stop("Missing runtime Stan compilation headers: ",
  paste(missing, collapse = ", "), ". Install BH and RcppEigen explicitly.")
stan_header_dirs <- c(boost_lib = rstan::rstan_options("boost_lib"),
  eigen_lib = rstan::rstan_options("eigen_lib"))
if (length(stan_header_dirs) != 2L || any(!nzchar(stan_header_dirs)) || any(!dir.exists(stan_header_dirs)))
  stop("RStan cannot locate its Boost/Eigen header directories; inspect its rstan_options().")
cat("Runtime Stan headers:\n")
print(headers)
cat("RStan header directories:\n")
print(stan_header_dirs)
