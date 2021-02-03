# Perform a *successful* step
function DiffEqBase.step!(integrator::ODEFilterIntegrator)
    loopheader!(integrator)
    check_error!(integrator) != :Success && return
    perform_step!(integrator)
    loopfooter!(integrator)
    while !integrator.accept_step
        loopheader!(integrator)
        check_error!(integrator) != :Success && return
        perform_step!(integrator)
        loopfooter!(integrator)
    end
end


function check_error!(integrator::ODEFilterIntegrator)
    code = check_error(integrator)
    if code != :Success
        integrator.retcode = code
        postamble!(integrator)
    end
    return code
end
