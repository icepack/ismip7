#!/bin/bash
# Shared environment for ISMIP7 icepack2 batch jobs, on any cluster.
# Sourced by the inversion and projection runners and by submit.sh. Nothing
# here submits anything.
#
# Two halves. site_core.sh answers the site questions (which cluster, where is
# Firedrake, what does the scheduler want, the per-user sites/local.env) and
# provides ismip7_activate, ismip7_chain_resources and ismip7_banner. This file
# adds the model configuration every inversion and projection shares. A job
# that brings its whole model configuration with it, as the timing campaign's
# lanes do, sources site_core.sh alone.
#
#   ISMIP7_SITE=iu_quartz sbatch ...     # name it
#   sbatch ...                           # or let the hostname choose
# shellcheck disable=SC1091
. "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/site_core.sh"

# --- repository and data ------------------------------------------------
# The gitignored artifacts a checkout does not carry -- the forcing tree and the
# observational rasters -- hang off the invoking checkout. A site whose copies
# live elsewhere names these two in its site file, which is sourced before these
# defaults apply.
export ISMIP7_DATA_ROOT="${ISMIP7_DATA_ROOT:-$ISMIP7_REPO/ISMIP7/AIS}"
# simulation.py and preflight.py default this to the running checkout's
# antarctica/data, which is gitignored and empty in a fresh clone, so state it
# here too: the shell layer and the Python layer then name the same rasters.
export ISMIP7_OBS_DATA_ROOT="${ISMIP7_OBS_DATA_ROOT:-$ISMIP7_REPO/antarctica/data}"

# --- model configuration shared by every run ----------------------------
# The production mesh is 1000 m / 10 km with a 20 km buffer, the submission
# mesh since 25 September 2026 (issue 20): the finest pair the timing matrix
# carries through a 285-year run in two days at dt 0.05
# (TIMING_MATRIX_QUARTZ_SCPC_GAMG.md, 10.0 min per simulated year on 64
# ranks), and about four at the production step, 0.025, which
# projection.sbatch sets. The mesh file is Rice's build, the one the MAP is
# inverted on (antarctica/README.md, "Starting state"). The mesh name follows
# mesh_naming.mesh_basename, so naming another pair names its mesh too. The
# 2500 m configuration this replaces had a mesh outside that rule and is
# ISMIP7_LC=2500 ISMIP7_LC_COARSE=64000
# ISMIP7_MESH=$ISMIP7_REPO/antarctica/mesh/antarctica_64000_2500.msh.
#
# No solver is chosen here. The inversion sources this file as well, its
# linear solve is the full mixed-Jacobian MUMPS by construction, and the solver
# it finds in the environment is stamped on the MAP it writes; the forward
# runner (projection.sbatch) names its own.
#
# Every inversion runs regularized Coulomb. Budd's shelf gate was a sign test
# on the roundoff residue of N, so every Budd MAP predating that fix has to be
# re-inverted; set ISMIP7_FRICTION=budd explicitly for those re-inversions.
export ISMIP7_GEOMETRY_SPACE="${ISMIP7_GEOMETRY_SPACE:-dg0}"
export ISMIP7_FRICTION="${ISMIP7_FRICTION:-regularized_coulomb}"
export ISMIP7_N_FLOW="${ISMIP7_N_FLOW:-3.0}"
export ISMIP7_LC="${ISMIP7_LC:-1000}"
export ISMIP7_LC_COARSE="${ISMIP7_LC_COARSE:-10000}"
# int(float(buffer)), as mesh_naming.buffer_tag writes it.
_ismip7_buffer="${ISMIP7_BUFFER_M:-20000}"
export ISMIP7_MESH="${ISMIP7_MESH:-$ISMIP7_REPO/antarctica/mesh/antarctica_${ISMIP7_LC_COARSE}_${ISMIP7_LC}_buffered${_ismip7_buffer%%.*}.msh}"
# The MAP this configuration writes and reads: named from the law, so switching
# ISMIP7_FRICTION switches the file and a Budd re-inversion cannot land on the
# RC MAP. One variable for both halves of the workflow - inversion.sbatch
# writes it, projection.sbatch loads it - because the runners' `logvelnet`
# name is not one a forward can derive for itself.
. "$_ISMIP7_BR_DIR/../ismip7_names.sh"
export ISMIP7_MAP_DEFAULT="${ISMIP7_MAP_DEFAULT:-$ISMIP7_REPO/antarctica/mesh/$(ismip7_map_basename "$ISMIP7_FRICTION" "$ISMIP7_LC")}"

ismip7_banner_model() {
    echo "    mesh    $ISMIP7_MESH"
    echo "    lc=$ISMIP7_LC lc_coarse=$ISMIP7_LC_COARSE geometry=$ISMIP7_GEOMETRY_SPACE" \
         "friction=$ISMIP7_FRICTION n=$ISMIP7_N_FLOW"
}
