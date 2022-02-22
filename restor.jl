using Images
# using IO
using LinearAlgebra
using Colors
using DelimitedFiles
using Statistics
using ImageView
pillar = load("./gpillars.png")
graffiti = load("./ggraf.png")
pillar = Gray.(pillar)
graffiti = Gray.(graffiti)
# # ### pillar = load("/home/dsc/Dropbox/School/CompPhys/project-2/gpillars.png")
### graffiti = load("/smallerGraffiti.png")
global xDim = size(pillar)[2]
global yDim = size(pillar)[1]

"""Find the locations of the relevant pixels of the  graffiti mask."""
function findBadPixels(grafPattern)
    badPixels = Tuple{Int64, Int64}[]
    for i= 1:size(graffiti)[1], j=1:size(graffiti)[2]
        if graffiti[i, j] == 0
            push!(badPixels, (i,j))
        end
    end
    return badPixels
end

"""Applies the grattifi masks """
function graffitiIT(img, grafPattern)
    badPixels = findBadPixels(graffiti)
    for (i,j) in badPixels
        img[i, j] = 0
    end
    return img
end



# """Finds the pixels at the boundaries of the graffiti."""
# function findNotBadPixels(img, grafPoints)
#     boundInd = Tuple{Int64, Int64}[]
#     # Eliminating consideration for pixels at the edges
#
#     if i == 1
#         if j ==
#     for (i, j) in grafPoints
#         boundaries = [(i + J, j + J) for I=-1:1, J=-1:1]
#         boundSum = 0
#         for (k, l) in boundaries
#             boundSum += img[k, l]
#         end
#         if boundSum == 9.0
#             push!(boundInd, (i,j))
#         end
#     end
#     return boundInd
# end

"""Mapping between s and ij space."""
function s_ij_map(inside)
    # Make s-space and i,j space transitions
    s_to_ij = Dict{Int32, Tuple{Int32, Int32}}()
    ij_to_s = Dict{Tuple{Int32, Int32}, Int32}()
    s = 1
    for (i, j) in inside
        s_to_ij[s] = (i, j)
        ij_to_s[(i, j)] = s
        s += 1
    end
    return s_to_ij, ij_to_s
end

""""Finds neighbors for a given inside point and filters out points that are out of bounds."""
function find_neighbors((i, j))

    tmpNeighbors = [(i + I, j + J) for I=-1:1, J=-1:1 if abs(I) != abs(J)]
    neighbors = copy(tmpNeighbors)

    for i = 1:length(tmpNeighbors)
        if tmpNeighbors[i][1] > yDim || tmpNeighbors[i][1] == 1
            deleteat!(neighbors, i)
        end
        if tmpNeighbors[i][2] > yDim || tmpNeighbors[i][2] == 1
            deleteat!(neighbors, i)
        end
    end
    return neighbors
end


"""Matrix method to infer values of missing points."""
function makeMatrices(badImg, badCoordinates)
    inside = badCoordinates
    outside = Tuple{Int32, Int32}[]

    A = zeros(Int32, size(inside)[1], size(inside)[1])
    B = zeros(Float64, size(inside)[1])
    # Check for neighbors to see which pixels are outside the graffiti
    for (i, j) in inside
        neighbors = find_neighbors((i, j))
        for (k, l) in neighbors
            if badImg[i ,j] != badImg[k, l]
                push!(outside, (k, l))
            end
        end
    end

    ij_to_s = s_ij_map(inside)[2]

    # Find neighbors for every pixel inside
    for (i, j) in inside
        neighbors = find_neighbors((i,j))
        # Check if said neighbor is located outside
        true_false_list = neighbors .∈ [outside]
        # If neighbor is outside, add it to the RHS
        if any(x -> x == true, true_false_list)
            outside_neighbors = [neighbors[i] for i in findall(x -> x == true, true_false_list)]
            for (I, J) in outside_neighbors
                B[ij_to_s[(i, j)]] -= badImg[I, J]
            end
        end
        # Add 1 to the matrix for every neighbor inside
         A[ij_to_s[(i, j)], ij_to_s[(i, j)]] = -1*(4  )#)- size(findall(x -> x == true, true_false_list))[1])
        for (I, J) in neighbors
            if (I, J) in inside
                A[ij_to_s[(i, j)], ij_to_s[(I, J)]] = 1
            end
        end
    end
    return A, B
end

"""Function to restore the image by solving the matrices"""
function restore(A, B, grafPoints, badImg)


    A_inv = inv(A)
    U = A_inv * B
    goodImg = copy(badImg)

    s_to_ij = s_ij_map(grafPoints)[1]

    for s=1:size(U)[1]
        i = s_to_ij[s][1]
        j = s_to_ij[s][2]
        goodImg[i, j] = U[s]./maximum(U)
    end


    return goodImg

end

"""Calculating the error based on (1) from the project instructions."""
function calcDiscrepancy(originalImg, changedImg, grafPoints)
    meanI = mean([grafPoints[i, j] for (i, j) in grafPoints])
    Σ_χ = 0
    Σ_σ = 0
    for (i, j) in grafPoints
        Σ_χ = Σ_χ + (changedImg[i,j] - originalImg[i, j])^2
        Σ_σ = Σ_σ + (originalImg[i,j] - meanI)^2
    end

    σ_2 = Σ_σ/(length(grafPoints) - 1)
    χ_2 = Σ_χ/(length(grafPoints) * σ_2)
    print("\n Discrepancy score is $(χ_2). \n The variance is $(σ_2).")

    return χ_2, σ_2
end

badImage = graffitiIT(pillar, graffiti)
grafPoints = findBadPixels(graffiti)

# bounds = findNotBadPixels(badImage, grafPoints)
A, B = makeMatrices(badImage, grafPoints)
restoredImg = restore(A, B, grafPoints, badImage)
χ_2_orig, σ_2_orig = calcDiscrepancy(pillar, badImage)
χ_2_res, σ_2_res = calcDiscrepancy(pillar, restoredImg)
# s_to_ij = s_ij_map(grafPoints)


# function loopUntilGood(grafPoints, badImage, originalImg, threshold = 0.01)
#
#     A, B = makeMatrices(badImage, grafPoints)
#     restoredImg = restore(A, B, grafPoints, badImage)
#     χ_2, σ_2 = calcDiscrepancy(originalImg, restoredImg)
#
#     for i=1:1
#         badImage = restoredImg
#         A, B = makeMatrices(badImage, grafPoints)
#         restoredImg = restore(A, B, grafPoints, badImage)
#         χ_2, σ_2 = calcDiscrepancy(originalImg, restoredImg)
#         restoredImg
#
#     end
#
#     return restoredImg
# end

# restoredImg = loopUntilGood(grafPoints, badImage, pillar)
restoredImg
# writedlm("A_file.txt", restoredImg )
