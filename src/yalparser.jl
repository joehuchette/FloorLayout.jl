function yalparser(target::String)

    fp = open(target, "r")
    lines = readlines(fp)

    block_x = Float64[]
    block_y = Float64[]
    block_area = Float64[]
    block_name = String[]
    arr = Vector{String}[]
    ids = Set{String}()
    W = 0.0
    H = 0.0

    it = 1
    while it <= length(lines)
        line = rstrip(chomp(lines[it]),collect(";"))
        if contains(line, "MODULE ") # start of module
            spl = split(line)
            push!(block_name, spl[2])
            it += 1
            line = rstrip(chomp(lines[it]),collect(";"))
            spl = split(line)
            typ = spl[2]
            it += 1
            line = rstrip(chomp(lines[it]),collect(";"))
            spl = split(line)
            nums = float(spl[3:end])
            x = abs(nums[1] - nums[5])
            y = abs(nums[2] - nums[6])
            push!(block_x, mean([nums[1],nums[5]]))
            push!(block_y, mean([nums[2],nums[6]]))
            push!(block_area, x*y)
            if typ == "GENERAL"
                # do nothing
            elseif typ == "PARENT"
                W = y
                H = x
                arr = Array(Vector{String}, length(block_name)-1)
                while !contains(line, "NETWORK")
                    it += 1
                    line = lines[it]
                end
                it += 1
                line = rstrip(chomp(lines[it]),collect(";"))
                while !contains(line, "ENDNETWORK")
                    spl = split(line)
                    idx = find(block_name .== spl[2])[1]
                    arr[idx] = spl[3:end]
                    for elem in spl[3:end]
                        push!(ids, elem)
                    end
                    it += 1
                    line = rstrip(chomp(lines[it]),collect(";"))
                end
            end
        end
        it += 1
    end

    # compute nets now
    nets = Vector{Int}[]
    for elem in ids
        this_net = Int[]
        for it in 1:length(arr)
            if elem in arr[it]
                push!(this_net, it)
            end
        end
        push!(nets, this_net)
    end

    return (W,H), block_area[1:end-1], nets

end

function load_yal(filesource::String, β::Real)
    (W,H), a, nets = yalparser(joinpath(Pkg.dir("FloorLayout"), "examples", "data", filesource*".yal"))
    N = length(a)

    b = β*ones(N)
    wub = min(sqrt(a.*b), W)
    hub = min(sqrt(a.*b), H)
    wlb = a ./ hub
    hlb = a ./ wub

    c = zeros(N,N)
    for net in nets
        k = length(net)
        for i in net, j in net
            if i < j
                c[i,j] += 1/(k-1)
                c[j,i] += 1/(k-1)
            end
        end
    end

    Problem(N, W, H, wlb, hlb, wub, hub, a, b, c)
end
