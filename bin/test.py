import plot_nearby


def test_mwafluxes():
    # Need to find the actual fluxes ot change the 0 to
    assert(plot_nearby.read_MWA_fluxes("/data/MWA/", "J001513", "GLEAM J001513-472706")) == 0
