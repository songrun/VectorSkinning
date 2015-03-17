# Skinning Cubic Bézier Splines

This code implements the cubic Bézier spline portion of our SIGGRAPH Asia 2014 paper
[Skinning Cubic Bézier Splines and Catmull-Clark Subdivision Surfaces](http://cs.gmu.edu/~ygingold/splineskin/)
by [Songrun Liu](http://cs.gmu.edu/~sliu11/), [Alec Jacobson](http://www.cs.columbia.edu/~jacobson/), and [Yotam Gingold](http://cs.gmu.edu/~ygingold/).
The implementation consists of a web GUI that communicates with a Python back end.

## Installation

The software depends on Python 2.7+ and several modules:

- numpy
- scipy
- autobahn
- twisted
- cffi
- cvxopt (optional; disabled by default)

(Note: On OS X, you should not install modules into the built-in system Python. Instead,
install a copy of Python with a package manager such as
[homebrew](http://brew.sh/) or [macports](http://www.macports.org/), or [fink](http://www.finkproject.org/).)

You can install these modules with your favorite package manager, or using pip:

    pip install numpy scipy
    pip install autobahn
    pip install cffi
    pip install cvxopt

(If you don't have `pip`, install it using your OS package manager or get it from http://www.pip-installer.org/en/latest/installing.html .)

The software depends on several compiled components.
Binaries are provided in the repository for OS X.
See below for compilation instructions (for other platforms, or when making changes).

## Running

Run
    python web-gui.py open
to run the computation back-end and (the `open` argument) launches
the web GUI in your default browser.
(You can run the web GUI independently, and then later run `python web-gui.py`
and click the "Reconnect" button in the GUI to attach to it.)


## Usage

Drag and drop an SVG file onto the web GUI to load the SVG.
There are many examples in `web-gui-tests` and in `results`.
Not every SVG can be loaded.
(A good one to start with is Alec Jacobson's clam: `web-gui-tests/clam-songrun.svg`.)

Click on the "Handles" pane to switch to the mode that allows you to
add and manipulate handles. In this mode, clicking inside the SVG view adds a handle;
dragging an existing handle modifies the transform associated with the handle.

Click on the "Control point" pane to inspect and change constraints on
interpolated control points.

Various global options are available in the GUI under "Toggles".
More global options are hidden inside `parameters.py`.
You will have to restart the server component (`web-gui.py`) to see their effect.
You can restart the server independently of the GUI. Every time you restart the server,
simply reconnect to it with the "Reconnect" button.
(Note: Some server-side options in `parameters.py` are "defaults", so they won't
have an effect when using the "Reconnect" button. For example, `kG1andAconstraints`
has an effect only when parsing an SVG initially, and `kArcLengthDefault` has
a button in the GUI under "Toggles" which overrides it.)


## Compiling from scratch

Although largely written in Python, there are several optional compiled components required.
Compiled versions are checked into the repository, but you can compile them from scratch.

1. (required for BBW weights) [Triangle](http://www.cs.cmu.edu/~quake/triangle.html). There should be a `triangle` directory containing a compiled binary of `triangle`.
2. (required for BBW weights) [libigl](https://github.com/libigl/libigl). There should be a directory or symbolic link to the root of `libigl` in the `src` directory.
3. (required for BBW weights) `bbw_wrapper` contains a compiled dynamic library. `cd` into the directory and run the appropriate one-liner located at the top of `bbw_wrapper/bbw.py`.
4. (optional for comparison with Schneider 1990 Graphics Gems) `FitCurves` contains a compiled dynamic library. `cd` into the directory and run the appropriate one-liner located at the top of `FitCurves/FitCurves.py`.
5. (optional for barycentric interpolation; disabled by default) `raytri` contains a compiled dynamic library. `cd` into the directory and run the appropriate one-liner located at the top of `raytri/raytri.py`.
