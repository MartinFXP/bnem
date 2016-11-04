exportVars <-
function(type = "ga") { # type: ga, loc, ex
  if ("ga" %in% type) {
    return(list("computeScoreNemT1", "simulateStatesRecursiveAdd", "simulateStatesRecursive", "reduceGraph", "getNemFit", "checkCNOlist", "checkNEMlist", "computeFc",  "computeSm", "relFit", "sizeFac", "method", "removeCycles", "dnf2adj", "verbose", "getHierarchy", "absorption", "checkMethod"))
  }
  if ("loc" %in% type) {
    return(list("computeScoreNemT1", "simulateStatesRecursiveAdd", "simulateStatesRecursive", "reduceGraph", "getNemFit", "checkCNOlist", "checkNEMlist", "computeFc",  "computeSm", "relFit", "sizeFac", "method", "absorptionII", "max.steps", "max.time", "removeCycles", "node", "dnf2adj", "verbose", "absorpII", "draw", "getHierarchy", "absorption", "bitStrings", "checkMethod"))
  }
  if ("ex" %in% type) {
    return(list("computeScoreNemT1", "simulateStatesRecursiveAdd", "simulateStatesRecursive", "reduceGraph", "getNemFit", "checkCNOlist", "checkNEMlist", "computeFc",  "computeSm", "relFit", "sizeFac", "method", "removeCycles", "dnf2adj", "verbose", "getHierarchy", "absorption", "checkMethod", "NEMlist"))
  }
}
