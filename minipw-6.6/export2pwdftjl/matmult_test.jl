nij = 4
nab = 3
ngm = 5

qgm = zeros(ComplexF64, ngm, nij)
aux = zeros(ComplexF64, ngm, nab)
deeaux = zeros(Float64, nij, nab)

qgm[:,:] .= 1.1 + 2.3im
qgm[2,2] = 7.0 + 0.0*im

aux[:,:] .= 5.0 + 3.0im
aux[2,3] = 9.7

using DelimitedFiles
res = reshape(readdlm("fort.100"), size(deeaux))