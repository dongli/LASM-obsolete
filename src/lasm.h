#ifndef __lasm__
#define __lasm__

#include "lasm_commons.h"
#include "Advection/Tracer/Tracer.h"
#include "Advection/Tracer/TracerManager.h"
#include "Advection/AdvectionManager.h"
#include "Advection/TestCase/AdvectionTestCase.h"

#if defined USE_CARTESIAN_DOMAIN
#include "Advection/TestCase/CartesianRotationTestCase.h"
#elif defined USE_SPHERE_DOMAIN
#include "Advection/TestCase/DeformationTestCase.h"
#include "Advection/TestCase/SolidRotationTestCase.h"
#include "Advection/TestCase/BarotropicTestCase.h"
#endif

#endif
