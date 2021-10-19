# !/bin/python
# Code for all the priortransformation functions to input to ultranest
# TODO: update the internalbremss and curve models with correct ranges 

def singhomobremss(cube):
    params = cube.copy()
    params[0] = 10 ** (cube[0] * 3 - 1)
    params[1] = cube[1] * 20 - 10
    params[2] = cube[2] * 0.5
    # params[3] = (10 ** (cube[3] * 2)) * 2
    return params


def singinhomobremss(cube):
    params = cube.copy()
    params[0] = 10 ** (cube[0] * 3 - 1)
    params[1] = cube[1] * 20 - 10
    params[2] = cube[2] * 2 - 1
    params[3] = cube[3] * 0.5
    # params[4] = ((10 ** (cube[4] * 2))) * 2
    return params


def internalbremss(cube):
    params = cube.copy()
    params[0] = 10 ** (cube[0] * 3 - 1)
    params[1] = cube[1] * 20 - 10
    params[2] = cube[2] * 0.5
    return params


def singSSA(cube):
    params = cube.copy()
    params[0] = 10 ** (cube[0] * 3 - 1)
    params[1] = cube[1] * 20 - 10
    params[2] = cube[2] * 0.5
    return params


def singhomobremsscurve(cube):
    params = cube.copy()
    params[0] = 10 ** (cube[0] * 3 - 1)
    params[1] = cube[1] * 20 - 10
    params[2] = cube[2] * 0.5
    params[3] = cube[3]*4 - 2
    return params 


def singinhomobremsscurve(cube):
    params = cube.copy()
    params[0] = 10 ** (cube[0] * 3 - 1)
    params[1] = cube[1] * 20 - 10
    params[2] = cube[2] * 2 - 1
    params[3] = cube[3] * 0.5
    params[3] = cube[3]*4 - 2
    return params 


def singinhomobremssbreakexp(cube):
    params = cube.copy()
    params[0] = 10 ** (cube[0] * 3 - 1)
    params[1] = cube[1] * 20 - 10
    params[2] = cube[2] * 2 - 1
    params[3] = cube[3] * 0.3
    params[4] = ((10 ** (cube[4] * 2))) * 2
    return params


def singhomobremssbreakexp(cube):
    params = cube.copy()
    params[0] = 10 ** (cube[0] * 3 - 1)
    params[1] = cube[1] * 20 - 10
    params[2] = cube[2] * 0.5
    params[3] = (10 ** (cube[3] * 2)) * 2
    return params


def singSSAbreakexp(cube):
    params = cube.copy()
    params[0] = 10 ** (cube[0] * 3 - 1)
    params[1] = cube[1] * 20 - 10
    params[2] = cube[2] * 0.5
    params[3] = (10 ** (cube[3] * 2)) * 2
    return params