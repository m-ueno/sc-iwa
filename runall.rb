programs = ["jacobi2", "gauss2", "cg"]
sizes = ["900", "2500", "3600"]

programs.each do |prog|
  sizes.each do |size|
    p [prog, size]
    system("./#{prog} #{size}")
  end
end

