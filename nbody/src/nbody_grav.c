/*
 * Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
 * Copyright (c) 2010-2011 Rensselaer Polytechnic Institute
 * Copyright (c) 2010-2011 Matthew Arsenault
 *
 * This file is part of Milkway@Home.
 *
 * This file is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include "nbody_priv.h"
#include "nbody_util.h"
#include "nbody_grav.h"
#include "milkyway_util.h"

#include "../../third_party/Expanse/src/bfe.h" 
#include "string.h"

#define EXTERNAL_POTENTIAL_BFE 3

#ifdef _OPENMP
  #include <omp.h>
#endif /* _OPENMP */

/**
 * @brief Initializes the BFE model from a file. Call this during simulation setup.
 * @param filename Path to the BFE coefficient file.
 */
void nbGravInitBFE(NBodyCtx* ctx) {
    if (ctx->bfeModel) {
        bfe_destroy(ctx->bfeModel);
    }
    if (!ctx->coeffs_file_path)
        return;
    ctx->bfeModel = bfe_create_from_file(ctx->coeffs_file_path);
    if (!ctx->bfeModel) {
        mw_fail("Failed to initialize BFE model from file: %s\n", ctx->coeffs_file_path);
    }
}

/**
 * @brief Cleans up and frees the BFE model. Call this during simulation teardown.
 */
void nbGravCleanupBFE(NBodyCtx* ctx) {
    if (ctx->bfeModel) {
        bfe_destroy(ctx->bfeModel);
        ctx->bfeModel = NULL;
    }
}

// (The original nbGravity and nbGravity_Exact functions remain unchanged)
// ... [omitted for brevity] ...
static inline mwvector nbGravity(const NBodyCtx* ctx, NBodyState* st, const Body* p)
{
    mwbool skipSelf = FALSE;
    mwvector pos0 = Pos(p);
    mwvector acc0 = ZERO_VECTOR;
    const NBodyNode* q = (const NBodyNode*) st->tree.root;
    while (q != NULL) {
        mwvector dr = mw_subv(Pos(q), pos0);
        real drSq = mw_sqrv(dr);
        if (isBody(q) || (drSq >= Rcrit2(q))) {
            if (mw_likely((const Body*) q != p)) {
                real drab, phii, mor3;
                drSq += ctx->eps2;
                drab = mw_sqrt(drSq);
                phii = Mass(q) / drab;
                mor3 = phii / drSq;
                acc0.x += mor3 * dr.x;
                acc0.y += mor3 * dr.y;
                acc0.z += mor3 * dr.z;
                if (ctx->useQuad && isCell(q)) {
                    real dr5inv, drQdr, phiQ;
                    mwvector Qdr;
                    Qdr.x = Quad(q).xx * dr.x + Quad(q).xy * dr.y + Quad(q).xz * dr.z;
                    Qdr.y = Quad(q).xy * dr.x + Quad(q).yy * dr.y + Quad(q).yz * dr.z;
                    Qdr.z = Quad(q).xz * dr.x + Quad(q).yz * dr.y + Quad(q).zz * dr.z;
                    drQdr = Qdr.x * dr.x + Qdr.y * dr.y + Qdr.z * dr.z;
                    dr5inv = 1.0 / (sqr(drSq) * drab);
                    phiQ = 2.5 * (dr5inv * drQdr) / drSq;
                    acc0.x += phiQ * dr.x;
                    acc0.y += phiQ * dr.y;
                    acc0.z += phiQ * dr.z;
                    acc0.x -= dr5inv * Qdr.x;
                    acc0.y -= dr5inv * Qdr.y;
                    acc0.z -= dr5inv * Qdr.z;
                }
            } else {
                skipSelf = TRUE;
            }
            q = Next(q);
        } else {
             q = More(q);
        }
    }
    if (!skipSelf) {
        nbReportTreeIncest(ctx, st);
    }
    return acc0;
}
static mwvector nbGravity_Exact(const NBodyCtx* ctx, NBodyState* st, const Body* p)
{
    int i;
    const int nbody = st->nbody;
    mwvector a = ZERO_VECTOR;
    const real eps2 = ctx->eps2;
    for (i = 0; i < nbody; ++i) {
        const Body* b = &st->bodytab[i];
        if (p == b) continue;
        mwvector dr = mw_subv(Pos(b), Pos(p));
        real drSq = mw_sqrv(dr) + eps2;
        real drab = mw_sqrt(drSq);
        real phii = Mass(b) / drab;
        real mor3 = phii / drSq;
        mw_incaddv(a, mw_mulvs(dr, mor3));
    }
    return a;
}

static mwvector bfe_grav(Body* b, BFEModel* bfe_model) {
    double bfe_pos[3];
    double bfe_force[3];
    double bfe_potential; // required by function, but not used here
    mwvector total_acc;

    mwvector p_vec = Pos(b);
    bfe_pos[0] = p_vec.x;
    bfe_pos[1] = p_vec.y;
    bfe_pos[2] = p_vec.z;

    // Call the external function, using the global BFE model pointer.
    bfe_calculate_potential_and_force(bfe_pos, bfe_model, &bfe_potential, bfe_force);
    
    SET_VECTOR(total_acc, bfe_force[0], bfe_force[1], bfe_force[2]);
    return total_acc;
}


static inline void nbMapForceBody(const NBodyCtx* ctx, NBodyState* st)
{
    int i;
    const int nbody = st->nbody;
    mwvector LMCx;
    mwvector a, externAcc;
    const Body* b;
    real lmcmass, lmcscale;

    const Body* bodies = mw_assume_aligned(st->bodytab, 16);
    mwvector* accels = mw_assume_aligned(st->acctab, 16);
    real barTime = st->step * ctx->timestep - st->previousForwardTime;

    if (ctx->LMC) {
        LMCx = st->LMCpos;
        lmcmass = ctx->LMCmass;
        lmcscale = ctx->LMCscale;
    } else {
        SET_VECTOR(LMCx,0.0,0.0,0.0);
        lmcmass = 0.0;
        lmcscale = 1.0;
    }

    if (ctx->potentialType == EXTERNAL_POTENTIAL_BFE && !((NBodyCtx*)ctx)->bfeModel)
        nbGravInitBFE((NBodyCtx*)ctx);

  #ifdef _OPENMP
    #pragma omp parallel for private(i, b, a, externAcc) shared(bodies, accels) schedule(dynamic, 4096 / sizeof(accels[0]))
  #endif
    for (i = 0; i < nbody; ++i)
    {
        switch (ctx->potentialType)
        {
            case EXTERNAL_POTENTIAL_DEFAULT:
                b = &bodies[i];
                a = nbGravity(ctx, st, b);
                externAcc = mw_addv(nbExtAcceleration(&ctx->pot, Pos(b), barTime), plummerAccel(Pos(b), LMCx, lmcmass, lmcscale));
                mw_incaddv(a, externAcc);
                accels[i] = a;
                break;
            case EXTERNAL_POTENTIAL_NONE:
                accels[i] = nbGravity(ctx, st, &bodies[i]);
                break;
            case EXTERNAL_POTENTIAL_CUSTOM_LUA:
                a = nbGravity(ctx, st, &bodies[i]);
                nbEvalPotentialClosure(st, Pos(&bodies[i]), &externAcc);
                mw_incaddv(externAcc, plummerAccel(Pos(&bodies[i]), LMCx, lmcmass, lmcscale));
                mw_incaddv(a, externAcc);
                accels[i] = a;
                break;
            case EXTERNAL_POTENTIAL_BFE:
                b = &bodies[i];
                accels[i] = bfe_grav(b, ctx->bfeModel);
                break;
            default:
                mw_fail("Bad external potential type: %d\n", ctx->potentialType);
        }
    }

    if (ctx->potentialType == EXTERNAL_POTENTIAL_BFE)
        nbGravCleanupBFE((NBodyCtx*)ctx);
}


static inline void nbMapForceBody_Exact(const NBodyCtx* ctx, NBodyState* st)
{
    int i;
    const int nbody = st->nbody;
    mwvector LMCx;
    mwvector a, externAcc;
    const Body* b;
    real lmcmass, lmcscale;

    Body* bodies = mw_assume_aligned(st->bodytab, 16);
    mwvector* accels = mw_assume_aligned(st->acctab, 16);
    real barTime = st->step * ctx->timestep - st->previousForwardTime;

    if (ctx->LMC) {
        LMCx = st->LMCpos;
        lmcmass = ctx->LMCmass;
        lmcscale = ctx->LMCscale;
    } else {
        SET_VECTOR(LMCx,0.0,0.0,0.0);
        lmcmass = 0.0;
        lmcscale = 1.0;
    }
    
    if (ctx->potentialType == EXTERNAL_POTENTIAL_BFE && !((NBodyCtx*)ctx)->bfeModel)
        nbGravInitBFE((NBodyCtx*)ctx);

  #ifdef _OPENMP
    #pragma omp parallel for private(i, b, a, externAcc) shared(bodies, accels) schedule(dynamic, 4096 / sizeof(accels[0]))
  #endif
    for (i = 0; i < nbody; ++i)
    {
        switch (ctx->potentialType)
        {
            case EXTERNAL_POTENTIAL_DEFAULT:
                b = &bodies[i];
                a = nbGravity_Exact(ctx, st, b);
                externAcc = mw_addv(nbExtAcceleration(&ctx->pot, Pos(b), barTime), plummerAccel(Pos(b), LMCx, lmcmass, lmcscale));
                mw_incaddv(a, externAcc);
                accels[i] = a;
                break;
            case EXTERNAL_POTENTIAL_NONE:
                accels[i] = nbGravity_Exact(ctx, st, &bodies[i]);
                break;
            case EXTERNAL_POTENTIAL_CUSTOM_LUA:
                a = nbGravity_Exact(ctx, st, &bodies[i]);
                nbEvalPotentialClosure(st, Pos(&bodies[i]), &externAcc);
                mw_incaddv(externAcc, plummerAccel(Pos(&bodies[i]), LMCx, lmcmass, lmcscale));
                mw_incaddv(a, externAcc);
                accels[i] = a;
                break;
            case EXTERNAL_POTENTIAL_BFE:
                b = &bodies[i];
                accels[i] = bfe_grav(b, ctx->bfeModel);
                break;
            default:
                mw_fail("Bad external potential type: %d\n", ctx->potentialType);
        }
    }

    if (ctx->potentialType == EXTERNAL_POTENTIAL_BFE)
        nbGravCleanupBFE((NBodyCtx*)ctx);
}

static inline NBodyStatus nbIncestStatusCheck(const NBodyCtx* ctx, const NBodyState* st)
{
    if (st->treeIncest) {
        return ctx->allowIncest ? NBODY_TREE_INCEST_NONFATAL : NBODY_TREE_INCEST_FATAL;
    }
    return NBODY_SUCCESS;
}

NBodyStatus nbGravMap(const NBodyCtx* ctx, NBodyState* st)
{
    if (ctx->potentialType == EXTERNAL_POTENTIAL_BFE) {
        if (!ctx->bfeModel) {
            mw_fail("BFE potential selected, but model is not initialized. Call nbGravInitBFE() first.\n");
        }
        // The 'Exact' mapping function is fine to use here, as it's just a loop.
        nbMapForceBody_Exact(ctx, st);

        if (st->potentialEvalError) {
            return NBODY_LUA_POTENTIAL_ERROR;
        }
        return NBODY_SUCCESS;
    }

    // Original logic for all other potential types
    NBodyStatus rc;
    if (mw_likely(ctx->criterion != Exact)) {
        rc = nbMakeTree(ctx, st);
        if (nbStatusIsFatal(rc))
            return rc;
        nbMapForceBody(ctx, st);
    } else {
        nbMapForceBody_Exact(ctx, st);
    }

    if (st->potentialEvalError) {
        return NBODY_LUA_POTENTIAL_ERROR;
    }
    return nbIncestStatusCheck(ctx, st);
}
