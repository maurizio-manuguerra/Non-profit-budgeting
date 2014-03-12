data.file=NULL #replace with a file name, e.g. data.file="budget.csv"
workdir="~/Dropbox/0.Research/budget/"
n.proj=30
n.comm=5
budget = 50000
n.frontier.pts=8
prop.frontier.pts = 0.05
benefit.at.risk=0.02
degree.pref.fun = 5
seed=6
########################################
########################################
generate.data=function(){
	if (!is.null(seed)) set.seed(seed)
	proj.names <<- paste("P",1:n.proj,sep="")
	m <<- floor(runif(n.proj,-10000,2000)) 	
	if (sum(-m) <= budget) {
		cat("You can fund all the projects!\n")
	}
	score=-(1:length(m))
	score.var  <<-  floor(runif(n.proj,0,6))
	cscore = matrix(nrow=n.proj,ncol=n.comm)
	for (i in 1:n.comm){
		cscore[,i]=rank.sample(score,score.var)
	}
	cscore <<- as.data.frame(cscore)
	V=cbind(n.proj-cscore+1,m)
	names(V)=c(paste("c",1:n.comm,sep=""),"Return [$]")
	row.names(V)=proj.names
	print(as.data.frame(V))
}
########################################
optimum.budget=function(){
	var_pref=array()
	for (i in 1:n.proj) var_pref[i] = var(predict(fit5, data.frame(score.frontier.proj=as.numeric(cscore[i,]))))
	err_pref <<- sqrt(var_pref)
	mission.risk.limit = (benefit.at.risk*sum(preferences[pre.approved.proj.ind]))^2
	res = opt_budget(m,preferences,var_pref,budget,mission.risk.limit)
	for (i in 1:10){
		mission.risk.limit = (benefit.at.risk*res$benefit)^2
		res=opt_budget(m,preferences,var_pref,budget,mission.risk.limit)
		res$uncert=sqrt(res$uncert)
	}
	return(res)
}
########################################
display.optim.results=function(){
	funded.proj=res$funded.proj
	unfunded.proj=res$unfunded.proj
	n.funded.proj=res$n.funded.proj
	expenditure=res$expenditure
	benefit=res$benefit
	uncert=res$uncert
	cat("--------------------------------------\n")
	cat(" \n")
	cat(" \n")
	cat(paste("Funded",n.funded.proj,"initiatives over", n.proj,"\n"))
	cat(" \n")
	cat(paste("Initiative","    ","Return","    ","   Preference\n",sep=''))
	for (i in funded.proj) cat(paste("    ",proj.names[i],"      ",m[i],"    ",preferences[i],"\n",sep=''))
	cat(" \n")
	cat(paste(" Total return:   ",-round(expenditure,0),"[$]\n",sep=''))
	cat(paste(" Total benefit:   ",-round(benefit,0),"[$]\n",sep=''))
	cat(paste(" Total uncertainty on benefit   ",round(uncert,0),"[$], i.e. ",round(100*uncert/benefit,2),"% of the total benefit\n",sep=''))
	myplot(m,score,pch=19,main="Funded initiatives are in red.",xlab="Return [$]",ylab="Rank  ")
	points(m[funded.proj],score[funded.proj],pch=19,col='red')
	ind=order(preferences)
	lines(preferences[ind],score[ind],col="blue")
	for (i in 1:n.proj) lines(c(m[i]-err_pref[i]/2,m[i]+err_pref[i]/2),c(score[i],score[i]))
}
########################################
rank.sample=function(real.rank,score.var){
	n=length(real.rank)
	probs=(n-(0:(n-1)))/sum(1:n) 
	probs=probs[-real.rank]
	proj.avail=1:n
	rank=NULL
	for (i in 1:n){
		for (bonus in 0:n){
			rr=-real.rank[i]
			vv=score.var[i]+bonus
			intvl=max(rr-vv,1):min(rr+vv,n)
			inds=which(proj.avail %in% intvl)
			if (length(inds)>0) break()
		}
		if (length(inds)==1){
			rank[i]=proj.avail[inds]
		} else {
			rank[i]=sample(proj.avail[inds],1,prob=probs[inds])
		}
		yes.inds=which(proj.avail != rank[i])
		proj.avail=proj.avail[yes.inds]
		probs=probs[yes.inds]
	}
	return(rank)
}
########################################
opt_budget=function(m,preferences,var_pref,budget,mission.risk.limit){
	obj <- -preferences
	max <- TRUE
	types <- rep("I",length(obj))
	bounds <- list(upper = list(ind = (1:length(obj)), val = rep(1,length(obj))))
	rhs <- c(budget,mission.risk.limit)
	mat <- matrix(c(-m,var_pref),byrow=T,nrow=length(rhs),ncol=length(preferences))
	dir <- c("<=","<=")
	int.prog.fit = Rsymphony_solve_LP(obj, mat, dir, rhs, types = types, max = max, bounds = bounds)
	funded.proj = which(int.prog.fit$solution==1)
	unfunded.proj = which(int.prog.fit$solution==0)
	n.funded.proj = length(funded.proj)
	expenditure = sum(-m[funded.proj])
	benefit = -sum(preferences[funded.proj])
	uncert = sum(var_pref[funded.proj]) 	
	return(list(funded.proj=funded.proj,unfunded.proj=unfunded.proj,n.funded.proj=n.funded.proj,expenditure=expenditure,benefit=benefit,uncert=uncert))
}
########################################
myplot=function(...,label=T){
	plot(..., yaxt="n")
	axis(side=2,at=(1:n.proj),label=as.character(n.proj:1))
	if (label){
		labels=as.vector(paste("  ", proj.names, "  ", sep=""))
		pointLabel(m, apply(cscore,1,median), labels = labels, cex=0.5)
	}
}
########################################
require(MonoPoly)
require(Rsymphony)
setwd(workdir)
if (is.null(data.file)) { 
	generate.data()
} else {
	cscore <<- read.csv(data.file, header=TRUE)
	n.proj = nrow(cscore)
	n.comm = ncol(cscore)
	proj.names = row.names(cscore)
	if (n.proj<2 | n.comm<2) {
		cat("There are not enough projects or committee members in the data set.\n")
		return()
	}
}
score = apply(cscore,1,median)	
max.x=max(m)
min.x=min(m)
max.y=max(cscore)
min.y=min(cscore)
center.x=(max.x+min.x)/2
center.y=(max.y+min.y)/2
slope=-(max.y-min.y)/(max.x-min.x)
nsteps=100
step.y=2*(max.y-min.y)/(nsteps)
score.y.max=max.y+(m-max.x)*slope
c0=0 
for (i in 1:nsteps){
	score.limit = score.y.max -i*step.y
	funded.proj = which(score>score.limit)
	n.funded.proj = length(funded.proj)
	if (n.funded.proj > 0) {
		cost.funded.proj = sum(m[funded.proj])
	} else {
		cost.funded.proj = 0
	}
	c0=cost.funded.proj		
	if (-cost.funded.proj > budget) {
		score.limit = score.limit + step.y/2 
		break
	}
}
pre.approved.proj.ind = funded.proj 
for (i in 1:nsteps){
	score.int.down = score.limit - i*step.y
	score.int.up   = score.limit + i*step.y
	frontier.proj  = which(score>score.int.down & score<score.int.up)
	n.frontier.proj = length(frontier.proj)
	if (n.frontier.proj > n.frontier.pts & n.frontier.proj/n.proj > prop.frontier.pts) break
}
ind=order(m[frontier.proj])
m.frontier.proj = m[frontier.proj][ind] 
score.frontier.proj = cscore[frontier.proj,][ind,]
score.frontier.proj <- unlist(score.frontier.proj)
m.frontier.proj <- rep(m.frontier.proj,n.comm)
ind=order(m.frontier.proj)
m.frontier.proj=m.frontier.proj[ind]
score.frontier.proj=score.frontier.proj[ind]
fit5 <- monpol(m.frontier.proj ~ score.frontier.proj, degree=degree.pref.fun)
preferences = as.vector(predict(fit5, data.frame(score.frontier.proj=score)))

benefit=NULL
risk=NULL
top.risk=seq(0.001,0.05,.001)
for (benefit.at.risk in top.risk){
	res=optimum.budget()
	benefit=c(benefit,res$benefit)
	uncert=res$uncert
	risk=c(risk,res$uncert/res$benefit)
}
plot(risk*100,benefit,xlab="Total variance of the total benefit to the university",ylab="Total contribution to the mission [$]", main="Choose the preferred level of uncertainty and re-run the analyses.")
lines(risk*100,benefit)
benefit.at.risk=as.real(readline(prompt="\nBased on the plot, which percentage on the total benefit to the organization do you want to accept as maximum uncertainty? Please insert a number between 0 and 100 (omit the '%' sign).\nUncertainty = "))/100
res=optimum.budget()
display.optim.results()


