function top(nelx, nely, volfrac, penalMax, rmin)
    # input: >> top(60,20,0.5,3,3)
    # Material properties
    E0 = 1
    Emin = 1e-9
    nu = 0.3
    penal = 0.96
    # Prepare finite element analysis
    A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12]
    A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6]
    B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4]
    B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2]
    KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11])
    nodenrs = reshape([1:(1+nelx)*(1+nely)],1+nely,1+nelx)
    edofVec = reshape(2*nodenrs[1:end-1,1:end-1]+1,nelx*nely,1)
    edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1)
    iK = int(reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1))
    jK = int(reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1))
    # DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
    F = sparse([2],[1],[-1.0],2*(nely+1)*(nelx+1),1)
    U = zeros(2*(nely+1)*(nelx+1),1)
    fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)])
    alldofs = [1:2*(nely+1)*(nelx+1)]
    freedofs = setdiff(alldofs,fixeddofs)
    # Prepare filter
    iH = int(ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1))
    jH = int(ones(size(iH)))
    sH = zeros(size(iH))
    k = 0;
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (i1-1)*nely+j1
            for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                    e2 = (i2-1)*nely+j2
                    k = k+1
                    iH[k] = e1
                    jH[k] = e2
                    sH[k] = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2))
                end
            end
        end
    end
    H = sparse(vec(iH),vec(jH),vec(sH))
    Hs = sum(H,2)
    # Initialize iteration
    x = volfrac*ones(nely,nelx)
    xPhys = copy(x)
    loop = 0
    change = 1
    # Start iteration
    while change > 0.01
        loop += 1
        penal = min(penalMax, penal+0.04)
        # FE-ANALYSIS
        sK = reshape(KE[:]*(Emin+xPhys[:]'.^penal*(E0-Emin)),64*nelx*nely,1)
        K = sparse(vec(iK),vec(jK),vec(sK)); K = (K+K')/2
        # @time KK = cholfact(K[freedofs,freedofs]); U[freedofs] = KK \ F[freedofs]
        @time U[freedofs] = K[freedofs,freedofs] \ F[freedofs]
        # Objective function and sensitivity analysis
        ce = reshape(sum((U[edofMat]*KE).*U[edofMat],2),nely,nelx)
        c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce))
        dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce
        dv = ones(nely,nelx)
        dc[:] = H*(dc[:]./Hs)
        dv[:] = H*(dv[:]./Hs)
        # OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
        l1 = 0; l2 = 1e9; move = 0.2; xnew = 0
        while (l2-l1)/(l1+l2) > 1e-3
            lmid = 0.5*(l2+l1)
            xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))))
            xPhys[:] = (H*xnew[:])./Hs
            if sum(xPhys[:]) > volfrac*nelx*nely
                l1 = lmid
            else
                l2 = lmid
            end
        end
        change = max(abs(xnew[:]-x[:]))
        x = xnew
        # print result
        @printf(" It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n",loop,c, mean(xPhys[:]),change)
    end
end
