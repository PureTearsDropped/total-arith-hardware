using Quaternions, Octonions
for line in eachline(ARGS[1])
    v = parse.(Float64, split(line))
    if length(v) == 8
        a = Quaternion(v[1:4]...); b = Quaternion(v[5:8]...); c = a*b
        println(join((c.s, c.v1, c.v2, c.v3), " "))
    else
        a = Octonion(v[1:8]...); b = Octonion(v[9:16]...); c = a*b
        println(join((c.s, c.v1, c.v2, c.v3, c.v4, c.v5, c.v6, c.v7), " "))
    end
end
