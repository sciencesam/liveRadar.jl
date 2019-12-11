using ImageFiltering

function CFAR2D(inputmat,numGuardCells,numTrainCells,Pfa)

    # INPUTS:
    # inputmat = matrix of input image pixels
    # numGuardCells = number of guard cells
    # numTrainCells = number of training cells
    # Pfa = probability of false alarm

    if numGuardCells<1
        return inputmat
    end
    if numTrainCells<1
        return inputmat
    end

    filterSize = 1+numGuardCells+numTrainCells
    kernel = ones(filterSize,filterSize)
    kernel[numTrainCells:filterSize-numTrainCells,numTrainCells:filterSize-numTrainCells] .= 0
    N = length(filter(kernel->!iszero(kernel),kernel)) # number of non-zero elements
    kernel = (1/N) .* kernel

    alpha = N*(Pfa^(-1/N)-1) # multiplicative factor for thresholding

    Pn = imfilter(inputmat,kernel)  # matrix of noise estimates for each point in the input matrix

    T = alpha .* Pn # CFAR detection threshold (multiplicative factor times the noise estimate)

    output = T .- inputmat

    return output

end
