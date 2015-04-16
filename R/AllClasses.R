setClass(Class="Trajectories", representation=representation(path="character",filenames="character",trajlist="list",avd="numeric"))

setClass(Class="TransTrajectories",representation=representation(tmethod="character",breakpoints="list",tavd="numeric",ttrajlist="list",tfilenames="character",difftraj="list"),contains="Trajectories")

setClass(Class="SegTrajectories",representation=representation(sparams="character",smatrix="data.frame"),contains=c("Trajectories","TransTrajectories"))

setClass(Class="SegSeriesTrajectories",representation=representation(ct="numeric",ssmatrix="data.frame",ssparams="character"),contains=c("Trajectories","TransTrajectories","SegTrajectories"))
