## Don't do "from parameters import *", or you won't be able to see or make programmatic
## changes to these values.

EngineType = {
	'YSApproach' : 0,
	'FourControls' : 1,
	'TwoEndpoints' : 2,
	'Jacobian' : 3 }
	
kG1andAconstraints = True
kClampOn = True
kArcLengthDefault = False
kVerbose = 2
# kTransformControls = False
kNoOverlays = False		
## Show comparison curves to Schneider 1990?
kComputeComparisonCurves = False
kEngineType = EngineType['FourControls']