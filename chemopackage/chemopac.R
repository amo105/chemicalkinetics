rowSums.cpp <- function(D){
	as.vector(.Call( "rowSumscpp", D, PACKAGE = "chemopac" ))
}

d2Mhat5.cpp <- function(RESIDUAL, iaaA, D, D2, R){
	.Call( "d2Mhat5cpp", RESIDUAL, iaaA, D, D2, R, PACKAGE = "chemopac" )
}

dMhat5.cpp <- function(RESIDUAL, iaaA, D, D2, R){
	.Call( "dMhat5cpp", RESIDUAL, iaaA, D, D2, R, PACKAGE = "chemopac" )
}

Mhat5.cpp <- function(RESIDUAL, iaaA, D2, R){
	.Call( "Mhat5cpp", RESIDUAL, iaaA, D2, R, PACKAGE = "chemopac" )
}

makeaa.cpp <- function(THETA, DESIGN, R){
	.Call( "makeaacpp", THETA, DESIGN, R, PACKAGE = "chemopac" )
}

diags.cpp <- function(X, Y){
	.Call( "diagscpp", X, Y, PACKAGE = "chemopac" )
}

solve.cpp <- function(X){
	.Call( "solvecpp", X, PACKAGE = "chemopac" )
}

odes5.cpp <- function(t, y, p){
	as.vector(.Call( "odes5cpp", t, y, p, PACKAGE = "chemopac" ))
}

makeA5.cpp <- function(D1, D2, D3, D4, D5, R){
	.Call( "makeA5cpp", D1, D2, D3, D4, D5, R, PACKAGE = "chemopac" )
}

kron2cpp <- function(aa0, aa1, yy){
	as.vector(.Call( "KRON2cpp", aa0, aa1, yy, PACKAGE = "chemopac" ))
}

kron3cpp <- function(aa0, aa1, aa2, yy){
	as.vector(.Call( "KRON3cpp", aa0, aa1, aa2, yy, PACKAGE = "chemopac" ))
}

disc.cpp <- function(dstar, psi, da, iaa){
	.Call( "disccpp", dstar, psi, da, iaa, PACKAGE = "chemopac" )
}








