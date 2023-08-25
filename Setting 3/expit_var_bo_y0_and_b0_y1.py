def expit_var_bo_y0_and_b0_y1(
        constants,
        b0_x,
        b0_y0s,
        b0_y1s,
        b1_y1_x,
        b1_y0,
        b2_y0,
        b2_y1,
        p0,
        p1,
        n,
        n1,
        r1,
        r2,
        PET, n_sim=100_000, dirname=None):
    dirname = dirname or f"./expit_var_bo_y0_and_b0_y1{int(p0*100)}"
    # intercepts = intercepts or [-1 for i in range(len(constants))]
    for c, b0_y0, b0_y1 in zip(constants, b0_y0s, b0_y1s):
        print(
            f"R -q -e 'source(\"expit_var_bo_y0_and_b0_y1.R\"); " + \
            f"sim(" + \
                f"{c}, {b0_x}, {b0_y0}, {b1_y1_x}, {b1_y0}, {b2_y0}, " + \
                f"{b2_y1}, {p0}, {p1}, {n}, {n1}, {r1}, {r2}, " + \
                f"{PET}, {b0_y1}, {n_sim}, \"{dirname}\")'"
        )

if __name__ == "__main__":
    expit_var_bo_y0_and_b0_y1(
        [0, 1, 2, 3, 4, 5],
        -1.18,  # b0_x
        [-1.1, -1.195, -1.315, -1.46, -1.625, -1.8], # b0_y0
        [-0.849, -1.18, -1.62, -2.1, -2.61, -3.11], # b0_y1
        .5,  # b1_y1 <- b1_x
        1/6,  # b1_y0
        1/3,  # b2_y0
        1, # b2_y1
        .25,  # p0
        .30,  # p1
        522,  # n
        223,  # n1
        57,  # r1
        146,  # r2
        .6112  # PET(p0)
    )

    expit_var_bo_y0_and_b0_y1(
        [0, 1, 1.5, 2, 3],
        -1,    # b0_x
        [-0.49, -1, -1.295, -1.615, -2.265],     # b0_y0
        [-0.122, -1, -1.28, -1.53, -2.05],    # b0_y1
        2,     #b1_y1 <- b1_x
        1,     #b1_y0
        1,     #b2_y0
        3,     #b2_y1
        .38,   #p0
        .47,   #p1
        212,   #n
        89,    #n1
        36,    #r1
        91,    #r2
        .7226 #PET(p0)
    )

    expit_var_bo_y0_and_b0_y1(
        [0, 1, 2, 3, 4, 5],
        -1,  # b0_x
        [-.89, -1.01, -1.145, -1.315, -1.49, -1.68], # b0_y0
        [-0.49, -1.015, -1.58, -2.13, -2.67, -3.22], # b0_y1
        .7,  # b1_y1 <- b1_x
        1/5,  # b1_y0
        1/3,  # b2_y0
        1.5,  # b2_y1
        .29,  # p0
        .38,  # p1
        195,  # n
        69,  # n1
        21,  # r1
        66,  # r2
        .6593,  # PET(p0)
    )
