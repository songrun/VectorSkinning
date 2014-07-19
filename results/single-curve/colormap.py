from PIL import Image
from numpy import *
import matplotlib.cm

def apply_colormap( arr, colormap_name, vmin = None, vmax = None, clip = False ):
    '''
    Given a 2D numpy array of scalar values,
    maps them according to the colormap named 'colormap_name'.
    Optional parameters vmin and vmax are the minimum and maximum values
    (default arr.min() and arr.max()).
    Values outside will be given the colormap's clipping colors, unless the
    optional 'clip' parameter is set to True.
    '''
    
    assert len( arr.shape ) == 2
    
    cmap = matplotlib.cm.ScalarMappable( cmap = colormap_name )
    result = cmap.to_rgba( arr.ravel(), bytes = True ).reshape( ( arr.shape[0], arr.shape[1], -1 ) )
    return result

def main():
    import os, sys
    def usage():
        print >> sys.stderr, "Usage:", sys.argv[0], "path/to/input colormap_name path/to/output"
        sys.exit(-1)
    
    argv = sys.argv[1:]
    if len( argv ) not in (2,3): usage()
    
    infile = argv.pop(0)
    colormap_name = argv.pop(0)
    if len( argv ) == 0:
        outfile = os.path.splitext( infile )[0] + '-' + colormap_name + '.png'
    else:
        outfile = argv.pop(0)
    
    assert len( argv ) == 0
    
    if os.path.exists( outfile ):
        print >> sys.stderr, "Output path exists, won't clobber:", outfile
        usage()
    
    data = array( Image.open( infile ).convert('L') )
    color_mapped = apply_colormap( data, colormap_name )
    Image.fromarray( color_mapped ).save( outfile )

if __name__ == '__main__': main()
