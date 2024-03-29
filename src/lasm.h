#ifndef __lasm__
#define __lasm__

#include "lasm_commons.h"
#include "Advection/Tracer/Tracer.h"
#include "Advection/Tracer/TracerManager.h"
#include "Advection/AdvectionManager.h"

#ifndef LASM_IN_ACTION
#include "Advection/TestCase/AdvectionTestCase.h"
#if defined LASM_CARTESIAN_DOMAIN
#include "Advection/TestCase/CartesianRotationTestCase.h"
#elif defined LASM_SPHERE_DOMAIN
#include "Advection/TestCase/DeformationTestCase.h"
#include "Advection/TestCase/SolidRotationTestCase.h"
#include "Advection/TestCase/BarotropicTestCase.h"
#include "Advection/TestCase/TerminatorChemistryTestCase.h"
#endif
#endif

#endif
