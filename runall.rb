programs = ["jacobi2", "gauss-seidel2", "cg"]
sizes = ["900", "2500", "3600"]

programs.each do |prog|
    sizes.each do |size|
        p [prog, size]
        1.times {
            system("./#{prog}.out #{size}")
        }
    end
end

