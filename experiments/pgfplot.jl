p = plot(rand(10));
filename = "asdf.tex";

function pgfsave(p, filename)
    filename = abspath(filename)
    savefig(p, filename)
    lines = readlines(filename)
    line = join(lines, "\n")
    @assert occursin(r", width=\{[0-9]+\.[0-9]*mm\}", line)
    @assert occursin(r", height=\{[0-9]+\.[0-9]*mm\}", line)
    new_line = replace(line, r", width=\{[0-9]+\.[0-9]*mm\}" => ", width=\\figwidth")
    new_line = replace(new_line, r", height=\{[0-9]+\.[0-9]*mm\}" => ", height=\\figheight")
    @assert !occursin(r", height=\{[0-9]+\.[0-9]*mm\}", new_line)
    @assert !occursin(r", width=\{[0-9]+\.[0-9]*mm\}", new_line)
    write(filename, new_line)
end
