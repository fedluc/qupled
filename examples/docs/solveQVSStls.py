import numpy as np
import qupled.quantum as qpq

theta_values = [1.000]
rs_values = [1.0]
step = 0.1
rs_min = 0.1
nl = 200

for theta in theta_values:
    for rs_max in rs_values:
        # Define a qVSStls object to solve the qVS-STLS scheme (first state point)
        print("Computing the qVS-STLS scheme, state point: theta = %5.3f, rs = %5.3f" % (theta, rs_min))
        qvsstls = qpq.QVSStls(rs_min, 
                  theta,
                  mixing = 0.5,
                  resolution = 0.1,
                  cutoff = 10,
                  matsubara = nl,
                  couplingResolution = 0.1,
                  degeneracyResolution = 0.1,
                  alpha = [-0.2, 0.4],
                  errorIntegrals = 1e-5,
                  iterations = 100,
                  errorAlpha = 1e-3,
                  threads = 8,
                  fixed=f"adr_fixed_theta{theta:.3f}_matsubara{nl:f}.bin")

        # Solve the QSTLS scheme
        qvsstls.compute()

        # Rest state points use precomputed FreeEnergyIntegrand
        rs_iter = np.arange(rs_min + step, rs_max + step, step)
        for rs in rs_iter:
            print("Computing the qVS-STLS scheme, state point: theta = %5.3f, rs = %5.3f" % (theta, rs))
            # Setup a new VSStls simulation and use the free energy
            qvsstls.inputs.coupling = rs
            qvsstls.inputs.alpha = [0.1, 0.5]
            rs_prev = rs - step
            qvsstls.setFreeEnergyIntegrand(f"rs{rs_prev:.3f}_theta{theta:.3f}_QVSSTLS.h5")
            qvsstls.compute()

# Plot the results
qvsstls.plot(["ssf"])

