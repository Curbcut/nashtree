library(nashtreeInference)
root <- tempfile('ntf-writer-test-'); dir.create(root)
rejects <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), 'error'))
}
path <- file.path(root, 'forest.ntf')
nodes <- rbind(c(0, 1, 2, 1, .5, 0), c(-1, -1, -1, 0, 0, 1), c(-1, -1, -1, 0, 0, 2))
ntf_write(path, 1, 4, 3, c(0,1), c(0,3), nodes, '{}')
h <- ntf_open(path)
x <- matrix(c(.5, .5+.Machine$double.eps, NA, -Inf, Inf), ncol=1)
stopifnot(identical(ntf_predict_fixed(h,x,3),matrix(c(1,2,1,1,2),5,1)))
stopifnot(identical(ntf_predict_mapped(h,x,rep(3,5)),c(1,2,1,1,2)))
ntf_close(h)
rejects(ntf_predict_fixed(h,x,3))
rejects(ntf_write(path,1,4,3,c(0,1),c(0,3),nodes,'{}'))
for (field in 1:4) {
  bad <- nodes; bad[1,field] <- NA_real_
  rejects(ntf_write(file.path(root,'bad'),1,4,3,c(0,1),c(0,3),bad,'{}'))
}
rejects(ntf_write(file.path(root,'bad'),1,4,c(3,3),c(0,1),c(0,3),nodes,'{}'))
rejects(ntf_write(file.path(root,'bad'),1,4,3,c(0,1),c(0,4),nodes,'{}'))
if (.Platform$OS.type != 'windows') {
  ntf_sync(path); ntf_sync(root)
  rejects(ntf_sync(file.path(root,'missing')))
} else rejects(ntf_sync(path))
stopifnot(identical(list.files(root),'forest.ntf'))
unlink(root,recursive=TRUE)
cat('native array writer and durability boundary passed\n')
