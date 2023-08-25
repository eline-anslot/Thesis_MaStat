def functions(p_short, months_short, n_sim, p_LTE_1, n, n1, r1, r2, PET, method="Optimal"):
    for p in p_short:
        for time in months_short:
            print(f"R -q -e 'source(\"sim_equi_nextra_time_short_var.R\"); design <- data.frame(method = c(\"{method}\"), n = c({n}), n1 = c({n1}), r1 = c({r1}), r2 = c({r2}), PET = c({PET})) ; sim({p}, {time}, {n_sim}, {p_LTE_1}, design)'")

if __name__ == "__main__":
    functions(
            [0.15, 0.25, 0.35, 0.45, 0.65, 1],
            [0, 2*4, 4*4, 6*4, 9*4, 11*4],
            100_000,
            .45,
            41,
            17,
            5,
            14,
            0.7653
        )

    functions(
            [0.15, 0.25, 0.35, 0.45, 0.65, 1],
            [0, 2*4, 4*4, 6*4, 9*4, 11*4],
            100_000,
            .35,
            149,
            56,
            15,
            45,
            0.6853
        )
    
    functions(
            [0.15, 0.25, 0.30, 0.45, 0.65, 1],
            [0, 2*4, 4*4, 6*4, 9*4, 11*4],
            100_000,
            .30,
            522,
            223,
            57,
            146,
            0.6112
        )


    functions(
            [0.15, 0.25, 0.28, 0.45, 0.65, 1],
            [0, 2*4, 4*4, 6*4, 9*4, 11*4],
            100_000,
            .28,
            1366,
            749,
            191,
            367,
            0.6423
        )