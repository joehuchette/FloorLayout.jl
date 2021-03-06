export get_problem_data

function get_problem_data(T::Symbol, b::Real)
    Lˣ, Lʸ, α, c = problem_specific(Val{T}())
    N = length(α)
    β = b*ones(N)
    ubˣ = min(sqrt(α.*β), Lˣ)
    ubʸ = min(sqrt(α.*β), Lʸ)
    lbˣ = α ./ ubˣ
    lbʸ = α ./ ubʸ

    @assert N == length(α) == size(c,1) == size(c,2) == length(β) ==
            length(lbˣ) == length(lbʸ) == length(ubˣ) == length(ubʸ)
    @assert issym(c)

    Problem(N, Lˣ, Lʸ, lbˣ, lbʸ, ubˣ, ubʸ, α, β, c)
end

function problem_specific{T}(::Val{T}) # fallback to the MCNC benchmarks
    T in [:hp, :apte, :xerox, :ami33, :ami49] || error("Unrecognized benchmark $T")
    (W,H), α, nets = yalparser(joinpath(Pkg.dir("FloorLayout"),"data",string(T)*".yal"))
    N = length(α)
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
    return W, H, α, c
end

function problem_specific(::Val{:small_test})
    N = 2
    α = [.25, .3]#, .35]
    c = rand(N,N)
    c += c'
    return 1, 1, α, c
end

function problem_specific(::Val{:hp4})
    (W,H), α, nets = yalparser(joinpath(Pkg.dir("FloorLayout"),"examples","data","hp.yal"))
    N = 4
    c = zeros(N,N)
    for net in nets
        k = length(net)
        for i in net, j in net
            (i > N || j > N) && continue
            if i < j
                c[i,j] += 1/(k-1)
                c[j,i] += 1/(k-1)
            end
        end
    end
    return 100W, 100H, α[1:N], c
end

function problem_specific(::Val{:hp5})
    (W,H), α, nets = yalparser(joinpath(Pkg.dir("FloorLayout"),"examples","data","hp.yal"))
    N = 5
    c = zeros(N,N)
    for net in nets
        k = length(net)
        for i in net, j in net
            (i > N || j > N) && continue
            if i < j
                c[i,j] += 1/(k-1)
                c[j,i] += 1/(k-1)
            end
        end
    end
    return 100W, 100H, α[1:N], c
end

function problem_specific(::Val{:hp7})
    (W,H), α, nets = yalparser(joinpath(Pkg.dir("FloorLayout"),"examples","data","hp.yal"))
    N = 7
    c = zeros(N,N)
    for net in nets
        k = length(net)
        for i in net, j in net
            (i > N || j > N) && continue
            if i < j
                c[i,j] += 1/(k-1)
                c[j,i] += 1/(k-1)
            end
        end
    end
    return W, H, α[1:N], c
end

function problem_specific(::Val{:Armour62_1})
    Lˣ, Lʸ = 22, 33 # pad a bit
    α = [27, 18 ,27, 18, 18, 18, 9, 9, 9, 24, 60, 42, 18, 24, 27, 75, 64, 41, 27, 45]
    #     1  2  3  4  5  6  7  8  9 10 11 12
    c = [   0   1.8   1.2     0    0    0     0     0    0 1.04  1.12     0    0  1.20    0     0     0     0     0     0
          1.8     0  0.96 24.45 0.78    0 13.95     0 1.20 1.35     0     0    0     0    0     0     0     0  6.90     0
          1.2  0.96     0     0    0 2.21     0     0 3.15 3.90     0     0    0 13.05    0     0     0     0 13.65     0
            0 24.45     0     0 1.08 5.70  7.50     0 2.34    0     0  1.40    0     0    0     0     0  1.50 15.75     0
            0  0.78     0  1.08    0    0  2.25  1.35    0 1.56     0     0    0     0 1.35     0     0     0     0     0
            0     0  2.21  5.70    0    0  6.15     0    0    0     0  0.45    0     0    0     0     0  1.05     0     0
            0 13.95     0  7.50 2.25 6.15     0 24.00    0 1.87     0     0    0  0.96    0     0     0  1.65     0  3.75
            0     0     0     0 1.35    0  24.0     0    0    0     0     0 0.60     0    0     0     0     0  7.50 33.45
            0  1.20  3.15  2.34    0    0     0     0    0    0     0     0    0  7.50    0     0  7.50     0     0     0
         1.04  1.35  3.90     0 1.56    0  1.87     0    0    0  0.36  12.0    0  18.6 1.92     0     0     0  5.25     0
         1.12     0     0     0    0    0     0     0    0 0.36     0  2.25    0  3.00 0.96  2.25     0     0     0     0
            0     0     0  1.40    0 0.45     0     0    0 12.0  2.25     0    0     0 1.65     0 15.00     0  8.40     0
            0     0     0     0    0    0     0  0.60    0    0     0     0    0  8.00 1.04  6.00     0     0     0     0
         1.20     0 13.05     0    0    0  0.96     0 7.50 18.6  3.00     0 8.00     0 9.75     0     0  0.90     0     0
            0     0     0     0 1.35    0     0     0    0 1.92  0.96  1.65 1.04  9.75    0     0  5.25     0     0     0
            0     0     0     0    0    0     0     0    0    0  2.25     0 6.00     0    0     0 12.00     0     0     0
            0     0     0     0    0    0     0     0 7.50    0     0  15.0    0     0 5.25 12.00     0     0  7.50     0
            0     0     0  1.50    0 1.05  1.65     0    0    0     0     0    0  0.90    0     0     0     0  4.65     0
            0  6.90 13.65 15.75    0    0     0  7.50    0 5.25     0  8.40    0     0    0     0  7.50  4.65     0     0
            0     0     0     0    0    0  3.75 33.45    0    0     0     0    0     0    0     0     0     0     0     0]
    return Lˣ, Lʸ, α, c
end

function problem_specific(::Val{:Armour62_2})
    Lˣ, Lʸ = 22, 33 # pad a bit
    #       A    B    C    D    E    F    G    H    J    K    L    M    N    P    R    S    T    U    V    W
    α = [  27,  18,  27,  18,  18,  18,   9,   9,   9,  24,  60,  42,  18,  24,  27,  75,  64,  41,  27,  45]
    c = [   0  120   80    0    0    0    0    0    0   40   80    0    0   80    0    0    0    0    0    0  #  A
            0    0   80 1630   30    0  930    0   80   90    0    0    0    0    0    0    0    0  460    0  #  B
            0    0    0    0    0  130    0    0  210  260    0    0    0  870    0    0    0    0  910    0  #  C
            0    0    0    0   60  380  500    0  130    0    0   70    0    0    0    0    0  100 1050    0  #  D
            0    0    0    0    0    0  150   90    0   60    0    0    0    0   90    0    0    0    0    0  #  E
            0    0    0    0    0    0  410    0    0    0    0   30    0    0    0    0    0   70    0    0  #  F
            0    0    0    0    0    0    0 1600    0  110    0    0    0   60    0    0    0  110    0  250  #  G
            0    0    0    0    0    0    0    0    0    0    0    0   40    0    0    0    0    0  500 2230  #  H
            0    0    0    0    0    0    0    0    0    0    0    0    0  500    0    0  500    0    0    0  #  J
            0    0    0    0    0    0    0    0    0    0   30  800    0 1240  160    0    0    0  350    0  #  K
            0    0    0    0    0    0    0    0    0    0    0  150    0  200   80 1500  350   90    0    0  #  L
            0    0    0    0    0    0    0    0    0    0    0    0    0    0  110    0 1000    0  360    0  #  M
            0    0    0    0    0    0    0    0    0    0    0    0    0  500   40  500    0   40    0    0  #  N
            0    0    0    0    0    0    0    0    0    0    0    0    0    0  650    0    0   60    0    0  #  P
            0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0  350    0    0    0  #  R
            0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0 1000    0    0    0  #  S
            0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0  500    0  #  T
            0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0  320    0  #  U
            0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0  #  V
            0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0] #  W
    c = c + c'

    return Lˣ, Lʸ, α, c
end

function problem_specific(::Val{:Bazaraa75_1})
    Lˣ, Lʸ = 10, 6
    α = [  9,  8, 10,  6,  4,  3,  3,  4,  2,   2,  1,  1,  2]
    c = [  0 288 180  54  72 180  27  72  36    0   0   9   0;
           0   0 240  54  72  24  48 160  16   64   8  16   0;
           0   0   0 120  80   0  60 120  60    0   0  30   0;
           0   0   0   0  72  18  18  48  24   48  12   0   0;
           0   0   0   0   0  12  12  64  16   16   4   8   0;
           0   0   0   0   0   0  18  24   6   12   3   3   0;
           0   0   0   0   0   0   0   0   6    6   3   6   0;
           0   0   0   0   0   0   0   0  16   16  16   4   0;
           0   0   0   0   0   0   0   0   0    4   4   2   0;
           0   0   0   0   0   0   0   0   0    0   2   2   0;
           0   0   0   0   0   0   0   0   0    0   0   2   0;
           0   0   0   0   0   0   0   0   0    0   0   0   0;
           0   0   0   0   0   0   0   0   0    0   0   0   0]
    c = c + c'

    return Lˣ, Lʸ, α, c
end

function problem_specific(::Val{:Bazaraa75_2})
    Lˣ, Lʸ = 9, 7

    α = [  9,  8,  9, 10,  6,  3,  3,  3,  2,   3,  2,  1,  1,  1]
    c = [  0  72 162  90 108  27   0   0  18   27  18   0   0   0;
           0   0  72  80   0  48   0  48  32    0  16   8   0   0;
           0   0   0  45  54  27  27  27   0   27   0   9  18   0;
           0   0   0   0  30   0  30  30  20   0   20  10  10   0;
           0   0   0   0   0  18   0  18  12   18  24   0   0   0;
           0   0   0   0   0   0   9   9   0    0   6   6   6   0;
           0   0   0   0   0   0   0   9  12    9   6   3   0   0;
           0   0   0   0   0   0   0   0   6    9   0   3   0   0;
           0   0   0   0   0   0   0   0   0    6   4   6   2   0;
           0   0   0   0   0   0   0   0   0    0   6   3   6   0;
           0   0   0   0   0   0   0   0   0    0   0   2   0   0;
           0   0   0   0   0   0   0   0   0    0   0   0   4   0;
           0   0   0   0   0   0   0   0   0    0   0   0   0   0;
           0   0   0   0   0   0   0   0   0    0   0   0   0   0]
    c = c + c'

    return Lˣ, Lʸ, α, c
end

function problem_specific(::Val{:Camp91})
    Lˣ, Lʸ = 25, 51
    α = [238, 112, 160, 80, 120, 80, 60, 85, 221, 119]
    N = length(α)
    c = zeros(N,N)
    c[ 1, 6] = 218
    c[ 2, 6] = 148
    c[ 2, 9] = 296
    c[ 3, 4] =  28
    c[ 3, 5] =  70
    c[ 4, 6] =  28
    c[ 4, 7] =  70
    c[ 4, 8] = 140
    c[ 5, 8] = 210
    c[ 7,10] =  28
    c[ 8,10] = 888
    c[ 9,10] =  59.2
    c = c + c'

    return Lˣ, Lʸ, α, c
end


function problem_specific(::Val{:Bozer91})
    Lˣ, Lʸ = 10, 20
    α = [12, 7, 6, 5, 7, 22, 22, 13, 7, 22, 9, 13, 4, 17, 25]
    N = length(α)
    c = zeros(N,N)
    c[ 1,15] =  240
    c[ 2, 1] =  240
    c[ 3, 4] = 1200
    c[ 4,10] = 1200
    c[ 5,14] =  600
    c[ 6, 8] =  480
    c[ 7, 8] =  480
    c[ 8,15] =  120
    c[ 9,10] =  600
    c[10,12] =  600
    c[11, 7] =  480
    c[12,15] =  600
    c[13, 7] =  480
    c[14,12] =  600
    c[15, 2] =   10
    c[15, 3] =   25
    c[15, 5] =   25
    c[15, 6] =   40
    c[15, 9] =   25
    c[15,11] =   40
    c[15,13] =   20
    c = c + c'

    return Lˣ, Lʸ, α, c
end

function problem_specific(::Val{:Bozer97_1})
    Lˣ, Lʸ = 12, 13
    α = [16, 16, 16, 36, 36, 9, 9, 9, 9]
    #    1 2 3 4 5 6 7 8 9
    c = [0 0 0 5 5 0 0 0 1
         0 0 0 3 3 0 0 0 1
         0 0 0 2 2 0 0 0 1
         0 0 0 0 0 4 4 0 0
         0 0 0 0 0 3 0 0 4
         0 0 0 0 0 0 0 0 2
         0 0 0 0 0 0 0 0 1
         0 0 0 0 0 0 0 0 0
         0 0 0 0 0 0 0 0 0]
    c = c + c'

    return Lˣ, Lʸ, α, c
end

function problem_specific(::Val{:Bozer97_2})
    Lˣ, Lʸ = 6, 8
    α = [1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 16, 16]
    #     1  2  3  4  5  6  7  8  9 10 11 12
    c = [ 0  2  0  0 10  0  0  0  9  0  0  0 #  1
          0  0  0  5  0  0  0  0  0  0  7  2 #  2
          0  0  0  0  0  2  9  0  0  0  0  0 #  3
          0  0  0  0  0  0  0  0  0  3  0  3 #  4
          0  0  0  0  0  4  0  1  0  0  0  5 #  5
          0  0  0  0  0  0  0  0  4  0  0  0 #  6
          0  0  0  0  0  0  0  0  0  0  0  0 #  7
          0  0  0  0  0  0  0  0  0  0  5  0 #  8
          0  0  0  0  0  0  0  0  0  0  0  3 #  9
          0  0  0  0  0  0  0  0  0  0  0  0 # 10
          0  0  0  0  0  0  0  0  0  0  0  1 # 11
          0  0  0  0  0  0  0  0  0  0  0  0]# 12
    c = c + c'

    return Lˣ, Lʸ, α, c
end

# α = rand(1:N, N)
# Lˣ = ceil(Int, rand(.8:0.01:1.2) * sqrt(sum(α)))
# Lʸ = ceil(Int, sum(α) / Lˣ)
# c = (tmp = sprand(N,N,.2); tmp += tmp'; for i in 1:N; tmp[i,i] = 0; end; full(tmp))
# @show α;
# @show Lˣ, Lʸ;
# @show c;

function problem_specific(::Val{:Rand9})
    Lˣ, Lʸ = 6, 6
    α = [7,4,4,1,1,1,1,8,8]
    c = [0.0 0.7986094862350401 0.8499233254733563 0.0 0.7157553198075688 0.0 0.7437930539668685 0.0 1.175010054208596
 0.7986094862350401 0.0 0.0 0.0 0.0 1.527404963754088 0.0 0.8993465217540104 0.9770788914543318
 0.8499233254733563 0.0 0.0 0.6928155492041943 0.0 0.327312670502387 0.0 0.01939277023677377 0.9952866480280969
 0.0 0.0 0.6928155492041943 0.0 0.11064546027043032 0.2141258056898243 0.998823070162278 0.06454412519826769 0.4958584789801965
 0.7157553198075688 0.0 0.0 0.11064546027043032 0.0 0.0 0.0 0.16233570993620972 0.0014645204713570337
 0.0 1.527404963754088 0.327312670502387 0.2141258056898243 0.0 0.0 0.0 0.0 0.0
 0.7437930539668685 0.0 0.0 0.998823070162278 0.0 0.0 0.0 0.0 0.0
 0.0 0.8993465217540104 0.01939277023677377 0.06454412519826769 0.16233570993620972 0.0 0.0 0.0 0.4233288909381965
 1.175010054208596 0.9770788914543318 0.9952866480280969 0.4958584789801965 0.0014645204713570337 0.0 0.0 0.4233288909381965 0.0]

    return Lˣ, Lʸ, α, c
end


function problem_specific(::Val{:Rand10})
    Lˣ, Lʸ = 8, 6
    α =  [8,4,1,7,3,7,4,9,3,1]
    #     1  2  3  4  5  6  7  8  9 10
    c =  c = [0.0 0.0 0.9049260691773204 0.0 0.0 0.0 0.0 0.9989655886238922 0.0 0.0
 0.0 0.0 1.5694241789894503 0.00354884895732277 0.3137620799856242 0.0 0.0 0.0 0.0 0.2173791447506559
 0.9049260691773204 1.5694241789894503 0.0 0.9382338874211305 0.0 0.0 0.0 0.6139693075117343 0.0 0.0
 0.0 0.00354884895732277 0.9382338874211305 0.0 0.0 0.0 0.8913110827285136 0.0 0.0 0.3414856679085707
 0.0 0.3137620799856242 0.0 0.0 0.0 0.0 0.0 0.0 0.6328948674240107 0.0
 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.7456193867033032 0.1612888230451428 0.0
 0.0 0.0 0.0 0.8913110827285136 0.0 0.0 0.0 0.1620216504491383 0.0 0.11427952217846271
 0.9989655886238922 0.0 0.6139693075117343 0.0 0.0 0.7456193867033032 0.1620216504491383 0.0 0.0 0.17837661199235666
 0.0 0.0 0.0 0.0 0.6328948674240107 0.1612888230451428 0.0 0.0 0.0 0.9729596485553473
 0.0 0.2173791447506559 0.0 0.3414856679085707 0.0 0.0 0.11427952217846271 0.17837661199235666 0.9729596485553473 0.0]
    return Lˣ, Lʸ, α, c
end

function problem_specific(::Val{:Rand11})
    Lˣ, Lʸ = 9, 8
    α = [8,3,2,8,9,9,5,5,1,10,10]
    c = [0.0 0.0 0.6681298698645621 0.5469619335589575 0.0 0.6656772900595389 0.5422023238674989 0.0 0.7837115852348486 0.0 0.6796644127670064
 0.0 0.0 0.0 0.0 0.9431420891034161 0.0 0.0 1.6639710742325924 0.8441121778543434 0.0 0.0
 0.6681298698645621 0.0 0.0 0.0 0.8257319756858643 0.5105926660160405 0.5592266444786895 0.0 0.4075160359832646 0.0 0.0
 0.5469619335589575 0.0 0.0 0.0 0.0 0.0 0.6821101875680364 0.0 0.0 0.676801127801794 0.1920125959201131
 0.0 0.9431420891034161 0.8257319756858643 0.0 0.0 0.0 0.0 0.0 0.2776509135014058 0.0 0.0
 0.6656772900595389 0.0 0.5105926660160405 0.0 0.0 0.0 0.35272294703294826 0.0 0.0 0.0 0.0
 0.5422023238674989 0.0 0.5592266444786895 0.6821101875680364 0.0 0.35272294703294826 0.0 0.0 0.0 0.0 0.026042138065528686
 0.0 1.6639710742325924 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
 0.7837115852348486 0.8441121778543434 0.4075160359832646 0.0 0.2776509135014058 0.0 0.0 0.0 0.0 0.3023441261401367 0.0
 0.0 0.0 0.0 0.676801127801794 0.0 0.0 0.0 0.0 0.3023441261401367 0.0 0.0
 0.6796644127670064 0.0 0.0 0.1920125959201131 0.0 0.0 0.026042138065528686 0.0 0.0 0.0 0.0]
    Lˣ, Lʸ, α, c
end

function problem_specific(::Val{:Rand12})
    Lˣ, Lʸ = 10, 8
    α = [9,8,3,10,8,8,8,3,1,2,8,4]
    c = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3461964587900239 0.18380626866158423 0.0 0.37615042725365067
 0.0 0.0 0.0 0.0 0.7630077316908308 0.0 0.3712064996520412 0.39425245684528076 0.15004194657037773 0.0 0.0 0.0
 0.0 0.0 0.0 0.0 0.7528412087794389 0.0 0.0 0.6581231637329426 0.0 0.0 0.3608842601476572 0.0
 0.0 0.0 0.0 0.0 0.0 0.45365915531690204 0.0 0.0 0.37276051649622444 0.0 0.9459258931718393 0.0
 0.0 0.7630077316908308 0.7528412087794389 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.38318117651430694
 0.0 0.0 0.0 0.45365915531690204 0.0 0.0 0.0 0.2402561854030516 0.0 0.6852890352623946 0.36193208333680826 0.0
 0.0 0.3712064996520412 0.0 0.0 0.0 0.0 0.0 0.0 0.009034050537933602 0.0 0.8895683276627007 0.8646601311994306
 0.0 0.39425245684528076 0.6581231637329426 0.0 0.0 0.2402561854030516 0.0 0.0 0.6149644252149276 0.9786535832072509 0.0 0.0
 0.3461964587900239 0.15004194657037773 0.0 0.37276051649622444 0.0 0.0 0.009034050537933602 0.6149644252149276 0.0 0.0 0.3045850350930539 0.0
 0.18380626866158423 0.0 0.0 0.0 0.0 0.6852890352623946 0.0 0.9786535832072509 0.0 0.0 0.0 0.0
 0.0 0.0 0.3608842601476572 0.9459258931718393 0.0 0.36193208333680826 0.8895683276627007 0.0 0.3045850350930539 0.0 0.0 0.2818960645332804
 0.37615042725365067 0.0 0.0 0.0 0.38318117651430694 0.0 0.8646601311994306 0.0 0.0 0.0 0.2818960645332804 0.0]
    Lˣ, Lʸ, α, c
end

function problem_specific(::Val{:Rand13})
    Lˣ, Lʸ = 12, 9
    α = [8,13,9,12,3,8,2,9,12,10,5,1,11]
    c = [0.0 0.0 0.414286710242429 0.0 0.0 0.0 0.0 0.5332972399929927 0.7234154459659385 0.21173685682394439 0.0 0.0 0.0
 0.0 0.0 0.96343563706096 0.0 0.17619287941756667 0.0 0.0 0.0 0.0 0.0 1.7021339452640891 0.0 0.7813668667684333
 0.414286710242429 0.96343563706096 0.0 0.0 0.2552188035621441 0.0 0.0 0.9001704075319319 0.0 1.1436923143730524 0.0 0.0 0.0
 0.0 0.0 0.0 0.0 0.0 0.19568561775658289 0.048898206368805486 0.0 0.04373639235309823 0.4614440155085817 0.0 0.0 0.0
 0.0 0.17619287941756667 0.2552188035621441 0.0 0.0 0.0 0.0 0.8066651557067757 0.0 0.0 0.0 0.0 0.1546833575219424
 0.0 0.0 0.0 0.19568561775658289 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.8600801439975654 0.0
 0.0 0.0 0.0 0.048898206368805486 0.0 0.0 0.0 0.5899797726504763 0.699413350811348 0.0 0.0 0.8198801301148517 0.0
 0.5332972399929927 0.0 0.9001704075319319 0.0 0.8066651557067757 0.0 0.5899797726504763 0.0 0.0 0.0 0.7072561052323294 0.0 0.8862009474224881
 0.7234154459659385 0.0 0.0 0.04373639235309823 0.0 0.0 0.699413350811348 0.0 0.0 0.0 0.1691967772636247 0.0 0.9937467799481123
 0.21173685682394439 0.0 1.1436923143730524 0.4614440155085817 0.0 0.0 0.0 0.0 0.0 0.0 0.3242846674152646 0.0 0.8821774461706708
 0.0 1.7021339452640891 0.0 0.0 0.0 0.0 0.0 0.7072561052323294 0.1691967772636247 0.3242846674152646 0.0 0.0 0.0
 0.0 0.0 0.0 0.0 0.0 0.8600801439975654 0.8198801301148517 0.0 0.0 0.0 0.0 0.0 0.0
 0.0 0.7813668667684333 0.0 0.0 0.1546833575219424 0.0 0.0 0.8862009474224881 0.9937467799481123 0.8821774461706708 0.0 0.0 0.0]
    Lˣ, Lʸ, α, c
end
