include("./sigma0.jl")
include("./sigmax.jl")
include("./sigmaz.jl")
using LinearAlgebra

sigma0=cal_sigma0()
sigmax=cal_sigmax()
sigmaz=cal_sigmaz()


"""
o 2x2 matrix, 
promote to 2^nx2^n matrix, acting on the l-th qubit
o=sigmax
l=1
n=2
O=promote(sigmax,1,2)
reduce(O,[1])
"""
function promote(o,l,n)
    id1=Matrix(I,2^(l-1),2^(l-1))
    id2=Matrix(I,2^(n-l),2^(n-l))
    kron(id1,o,id2)
end

import Pkg
Pkg.add("Debugger")
using Debugger
    
"""
compute local reduced
idxs, a list a index that will be maintained
idx=1
idxs=[1]
O=promote(sigmax,1,2)
reduce(O,[1])
@enter reduce(O,[2])
idxs=[2]
o1=rand(2,2)
o2=rand(2,2)
O=kron(o1,o2)
o1check=reduce(O,[1])/tr(o2)
o1check-o1
o2check=reduce(O,[2])/tr(o1)
o2check-o2
"""
function reduce(O,idxs)
    n_idxs=length(idxs)
    n_total=trunc(Int,log(2,size(O)[1]))
    idxs_remain=collect(setdiff(Set(1:n_total),Set(idxs)))
    n_remain=n_total-n_idxs
    o=zeros(2^n_idxs,2^n_idxs)
    #i1=1; i2=2
    for i1 in 1:2^n_idxs
        for i2 in 1:2^n_idxs
            i1bin=idxTobin(i1,n_idxs)
            i2bin=idxTobin(i2,n_idxs)
            o[i1,i2]=0.0
            #j=1
            for j in 1:2^n_remain
                jbin=idxTobin(j,n_remain)
                i1fullidx=fullidx(n_total,i1bin,idxs,jbin,idxs_remain)
                i2fullidx=fullidx(n_total,i2bin,idxs,jbin,idxs_remain)
                o[i1,i2]+=O[binToidx(i1fullidx),binToidx(i2fullidx)]
            end            
        end
    end
    o
end

function fullidx(n_total,i1bin,idxs,jbin,idxs_remain)
    i1fulldx=Vector{Int}(undef,n_total)
    i1fulldx[idxs].=i1bin
    i1fulldx[idxs_remain].=jbin
    i1fulldx
end

"""
idx start from 1!!
idx, 1,..,2^n
idx=1
n=4
idxTobin(1,4)
idxTobin(2,4)
idxTobin(3,4)
binToidx(idxTobin(4,4))
sum(abs.([binToidx(idxTobin(i,4))-i for i in 1:2^4]))
"""

function idxTobin(idx,n)
    reverse(digits(idx-1,base=2,pad=n))
end

function binToidx(arr,s=0)
    v = 1
    for i in view(arr,length(arr):-1:1)
        s += v*i
        v <<= 1
    end 
    s+1
end


sigmaz1=promote(sigmaz,1,3)
sigmaz2=promote(sigmaz,2,3)
sigmaz3=promote(sigmaz,3,3)
sigmax1=promote(sigmax,1,3)
sigmax2=promote(sigmax,2,3)
sigmax3=promote(sigmax,3,3)
zz=sigmaz1*sigmaz2+sigmaz2*sigmaz3+sigmaz3*sigmaz1
x=sigmax1+sigmax2+sigmax3
j=1.0
h=0.01
ham=j*zz-h*x

val,vec=eigen(ham)

# ham*vec[:,1]-val[1]*vec[:,1]
rhoeff=vec[:,1:1]*transpose(vec[:,1:1])

reduce(rhoeff,[1,2])
reduce(rhoeff,[1])

