

meth_names = data.frame( old = c("make.Z"),
                         new = c("elec.data" ) )

#meth_names = mutate( meth_names,
#                     new = gsub( "_", ".", old, fixed = TRUE) )
meth_names

library( tidyverse)
meth_names$cmd = map2_chr( meth_names$old,
                           meth_names$new,
                           sprintf, fmt= "perl -pi.bak -e 's/%s/%s/g' *.R" )
meth_names

getwd()
#setwd( "../../R")
setwd( here::here( "R" ) )
map( meth_names$cmd, system )
setwd( here::here( "tests/testthat/" ) )
map( meth_names$cmd, system )


setwd( here::here( "tests/testthat/" ) )
