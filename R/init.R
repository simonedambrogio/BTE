.onAttach<-function(libname, pkgname){
  packageStartupMessage("Hello! This is the BTE package, verision ",
                        utils::packageDescription("BTE")$Version, ".")

 # rstan_options(auto_write = TRUE)
 # options(mc.cores = parallel::detectCores())
}
