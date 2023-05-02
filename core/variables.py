from .reco import get_ene

def get_energy_n_nonpion_pixels(pixels, n_pixels=5, verbose = False):    
    if len(pixels) > n_pixels:
        pixels = pixels[0: n_pixels]
        return sum(get_ene(pix) for pix in pixels)
    return 0.

