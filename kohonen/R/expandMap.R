expandMap <- function(kohobj) {
  ## define new maps
  whatmaps <- kohobj$whatmap
  oldxdim <- kohobj$grid$xdim
  oldydim <- kohobj$grid$ydim
  gr <- somgrid(xdim = oldxdim * 2,
                ydim = oldydim * 2,
                topo = kohobj$grid$topo,
                neighbourhood.fct = as.character(kohobj$grid$neighbourhood.fct),
                toroidal = kohobj$grid$toroidal)
  nhbrdist <- unit.distances(gr)

  ## put old codes in the new map
  ncodes <- gr$xdim * gr$ydim
  noldcodes <- oldxdim * oldydim
  colnrs <- (1:noldcodes) %% oldxdim
  colnrs[colnrs == 0] <- oldxdim ## start from 1
  rownrs <- floor(((1:noldcodes)-1) / oldxdim) ## start from 0
  newindices <- rownrs*4*oldxdim + colnrs * 2 - 1

  nvars <- sapply(kohobj$codes[whatmaps], ncol)
  codes <- lapply(seq(along = whatmaps),
                  function(ii) matrix(0, ncodes, nvars[ii]))
  for (j in seq(along = whatmaps)) 
    codes[[ j ]][newindices,] <- kohobj$codes[[ whatmaps[j] ]]
  
  ## for all empty units interpolate from neighbours
  emptyUnits <- (1:ncodes)[-newindices]
  for (i in emptyUnits) {
    closeones <- which(abs(nhbrdist[i, newindices] - 1) < .001)
    for (j in seq(along = whatmaps))
      codes[[j]][i,] <- colMeans(codes[[j]][newindices[closeones],,drop=FALSE])
  }
  
  kohobj$codes[whatmaps] <- codes
  kohobj$changes <- NULL
  kohobj$nhbrdist <- nhbrdist
  kohobj$grid <- gr
  
  kohobj
}
