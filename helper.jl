module Helper

    export generatePosDefMatrix, generateCompleatableMatrix


    function generatePosDefMatrix(n::Int64,rng)
        X = rand(rng,n,n)
        #X = 1/2*(X+X')
        #X = X + n*eye(n)
        X = X*X'
        return X
    end

    """
    generateCompleatableMatrix(n::Int64,density::Float64,rng)

    Generates a n x n-matrix X with specified density that is guaranteed to be completable to a
    positive definite matrix by selecting other values for the zeros. The matrix is created by starting
    with a positive definite matrix Y and then randomly overwriting off-diagonal elements with zeros until
    the desired density is achieved.
    """
    function generateCompleatableMatrix(n::Int64,density::Float64,rng)
        if density <= 0 || density > 1
            error("Density value has to be between 0 and 1")
        end

        # create a random dense positive definite matrix of dimension n
        X = generatePosDefMatrix(n,rng)
        numZeros = 0


        #determine roughly how many zeros the matrix should have
        desiredNumZeros = Int(floor((n^2-n)*(1-density)))
        # make number even
        desiredNumZeros = desiredNumZeros - mod(desiredNumZeros,2)


        # determine max size of the zero blocks
        maxBlockWidth = Int(floor(n/5))
        if maxBlockWidth == 0
            maxBlockWidth = 1
        end
        maxBlockSize = maxBlockWidth^2

        while numZeros < desiredNumZeros
            numZerosAvailable = (desiredNumZeros - numZeros) / 2
            if numZerosAvailable <= maxBlockSize
                maxSize = Int(floor(sqrt(numZerosAvailable)))
            else
                maxSize = maxBlockWidth
            end
            #randomly pick a blockwidth to delete
            bw = rand(1:maxSize)

            i_right = rand(1:n-1)
            j_top = rand(i_right+1:n)

            if bw == 1
                i_left = i_right
                j_btm = j_top
            else
                i_left = i_right-(bw-1)
                j_btm = j_top+(bw-1)

                if i_left < 1
                    i_left = 1
                end
                if j_btm > n
                    j_btm = n
                end

            end

            # check how may zeros are alreay in the considered block
            zerosPresent = size(find(x->x==0,X[i_left:i_right,j_top:j_btm]),1) * 2

            # overwrite block and transposed block
            X[i_left:i_right,j_top:j_btm] = 0
            X[j_top:j_btm,i_left:i_right] = 0

            # count number of nonzeros
            numZeros += (2*(i_right-i_left+1)*(j_btm-j_top+1) - zerosPresent)
        end
        return X

    end


end #Module